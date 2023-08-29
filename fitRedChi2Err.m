function [xFit, redChi2Fit, xErr, chiUpper, chiTrial, paramTrial, interpPts, slopes, intercepts, paramLower, paramUpper] = fitRedChi2Err(data, dataErr, model, x0, errPts, fact, offset)
%fitRedChi2Err Calculate errorbars on parameters fitted with a reduced
%chi-squared
%   Give data, errorbars for the data, a fitting model, initial fitting
%   values, the number of sampled error points, a factor, and an offset for
%   trial values. Return the fitted parameters, the reduced chi-squared,
%   the errorbars on each fitted parameter, the threshold gof, the gof from
%   each fit, the (iterated) parameter value from each fit, regression
%   slopes, regression intercepts, lower value,  and upper value. The
%   errorbar is determined by taking four points about the threshold gof
%   and getting the two intercepts with that threshold then subtracting the
%   x-axis values (of each iterated parameter) and dividing by 2. Note
%   that x0, factor, offset are rows. Don't use factors over 1.

redChi2 = @(fnObs, fnErr, fnCal, nParam) sum((fnObs-fnCal).^2./fnErr.^2, 'all')/(length(fnObs)-nParam);
obj = @(x) redChi2(data, dataErr, model(x), length(x0));
[xFit, redChi2Fit] = fminsearch(obj, x0');

% Use the fitted parameters above to generate trial parameters for the
% errorbar calculation with the input factors and terms in trialParam.
chiUpper = redChi2Fit*(1+1/(length(data)-length(xFit))); % Upper bound of the gof for determining the errorbar. Check this since dof changes and may be undefined.

paramTrial = nan(errPts, length(xFit));
for i = 1:length(xFit)
    paramTrial(:, i) = linspace(xFit(i)*(1-fact(i))-offset(i), xFit(i)*(1+fact(i))+offset(i), errPts);
end
xTrialFit = nan(errPts, length(xFit)-1, length(xFit)); % Each column is the parameter value after fitting, depth is going from one iterated parameter to the next while holding other constant.
chiTrial = nan(errPts, length(xFit)); % Chi2 for every fit for the errorbar
xErr = nan(size(xFit'));
slopes = nan(2, length(xFit)); % First element is for the lower-parameter, second upper
intercepts = nan(2, length(xFit));
interpPts = nan(4, length(xFit));
paramLower = nan(1, length(xFit));
paramUpper = nan(1, length(xFit));

for j = 1:length(xFit)
    initParam = xFit';  % Fitting to one fewer parameter, only send the initial values for the fitted parameters. Can make this more efficient.
    initParam(j) = [];
    if j==1
        modelCalcErr = @(x, iterParam) model([iterParam, x(j:length(xFit)-1)]); % Notation is fun(var, x(1), x(2), x(3), x(4)), fun(x(1), var, x(2), x(3), x(4)), etc.
    elseif j<length(xFit)
        modelCalcErr = @(x, iterParam) model([x(1:j-1), iterParam, x(j:length(xFit)-1)]);
    else
        modelCalcErr = @(x, iterParam) model([x(1:j-1), iterParam]);
    end
    
    % Fit to the peak
    for k = 1:errPts
        modelInputErr = @(x) modelCalcErr(x, paramTrial(k, j));
        obj = @(x) redChi2(data, dataErr, modelInputErr(x), length(x0)); % Using nParam=length(x0) rather than length(initParam) so gof's match when iterating parameter
        [xTrialFit(k,:,j), chiTrial(k,j)] = fminsearch(obj, initParam);
    end
    int1 = find(chiTrial(:, j)>chiUpper & paramTrial(:, j)<xFit(j), 1, 'last');
    int2 = find(chiTrial(:, j)<chiUpper & paramTrial(:, j)<xFit(j), 1, 'first');
    int3 = find(chiTrial(:, j)<chiUpper & paramTrial(:, j)>xFit(j), 1, 'last');
    int4 = find(chiTrial(:, j)>chiUpper & paramTrial(:, j)>xFit(j), 1, 'first');
    if isempty(int1) || isempty(int2) || isempty(int3) || isempty(int4)
        disp('Error too large.') % Could not find the intercepts required to calculate an error bar, likely because the parameter wasn't changed enough to pass the gof threshold.
    else
        interpPts(:, j) = [int1; int2; int3; int4]; % Largest parameter with unsuitable chi that's less than the optimal parameter, smallest parameter with suitable chi that's less than the optimal parameter, largest parameter with suitable chi that's greater than the optimal parameter, smallest parameter with unsuitable chi that's greater than the optimal parameter
        slopes(1, j) = (chiTrial(interpPts(2, j), j)-chiTrial(interpPts(1, j), j))/(paramTrial(interpPts(2, j), j)-paramTrial(interpPts(1, j), j));
        slopes(2, j) = (chiTrial(interpPts(4, j), j)-chiTrial(interpPts(3, j), j))/(paramTrial(interpPts(4, j), j)-paramTrial(interpPts(3, j), j));
        intercepts(1, j) = chiTrial(interpPts(1, j), j)-slopes(1, j)*paramTrial(interpPts(1, j), j);
        intercepts(2, j) = chiTrial(interpPts(3, j), j)-slopes(2, j)*paramTrial(interpPts(3, j), j);
        paramLower(j) = (chiUpper-intercepts(1, j))/slopes(1, j); % Greatest integrated intensity with unsuitable chi that's smaller than the optimal intensity
        paramUpper(j) = (chiUpper-intercepts(2, j))/slopes(2, j); % Smallest integrated intensity with unsuitable chi that's greater than the optimal intensity
        xErr(j) = abs(paramUpper(j)-paramLower(j))/2; % Size of the errorbar for integrated intensity
    end
end
end