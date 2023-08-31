function [specHeatCalc, specHeatMag, entropy, entropyErr, specHeat, specHeatErr, temperature, fval, chiUpper] = greaseCorrEntropy(specHeatTotalPre, specHeatTotalErrPre, specHeatAddendaPre, specHeatAddendaErrPre, temperaturePre, n, fitLimLow, fitLimHigh, debTemp0, alpha0, errPts, fact, offset)
%GREASECORRENTROPY Calculate magnetic entropy from heat capacity data
%   Given the measured heat capacity in J/mol-K, fit to a Debye model with
%   scale factor before subtracting grease in case some transferred to
%   sample. Return the magnetic contribution to the heat capacity and the
%   entropy in J/mol-K. n is the number of atoms per formula unit. Include
%   the fitting range, initial guesses for the Debye temperature and alpha.
%   ErrPts, fact, and offset determine the sampled range of values when
%   finding the uncertainties using fitRedChi2Err. Return sorted data.

N = 6.022e23; % Avogadro's number
kB = 1.380649e-23; % Boltzmann constant

% Sort data from low-T to high-T
[temperature, tIndPre] = sort(temperaturePre);
specHeatTotal = specHeatTotalPre(tIndPre);
specHeatTotalErr = specHeatTotalErrPre(tIndPre);
specHeatAddenda = specHeatAddendaPre(tIndPre);
specHeatAddendaErr = specHeatAddendaErrPre(tIndPre);

% Pick out data for the fitting
temperatureFit = temperature(temperature<=fitLimHigh & temperature>=fitLimLow);
specHeatTotalFit = specHeatTotal(temperature<=fitLimHigh & temperature>=fitLimLow);
specHeatTotalErrFit = specHeatTotalErr(temperature<=fitLimHigh & temperature>=fitLimLow);

% Fit for Debye temperature, obtaining the errorbars.
grease = @(T, alpha) (1+alpha).*interp1(temperature, specHeatAddenda, T); % Already interpolated to measurement by multiVu, but must be same size as Debye contribution when plugging in other temperatures
deb = @(T, TDeb) 9.*n.*N.*kB.*(T./TDeb).^3.*integralDeb(TDeb, T);
specHeatFn = @(T, TDeb, alpha) deb(T, TDeb)+grease(T, alpha); % Debye model plus grease background after applying scale factor
specHeatFnInput = @(x) specHeatFn(temperatureFit, x(1), x(2));
x0 = [debTemp0, alpha0];
[x, fval, xErr, chiUpper, chiTrial, paramTrial, interpPts, slopes, intercepts, paramLower, paramUpper] = fitRedChi2Err(specHeatTotalFit, specHeatTotalErrFit, specHeatFnInput, x0, errPts, fact, offset);

% Calculate the lattice contribution to specific heat at higher
% temperatures
specHeat = specHeatTotal-grease(temperature, x(2)); % Sample specific heat after subtracting grease from total
specHeatErr = sqrt(specHeatTotalErr.^2+specHeatAddendaErr.^2);
specHeatCalc = deb(temperature, x(1)); % Calculated sample specific heat

% Calculate the magnetic heat capacity
specHeatMag = specHeat-specHeatCalc;

% Calculate the magnetic entropy. See 04/06/2023 notes. The error in the
% magnetic component of the specific heat equals the error in the specific
% heat.
entropy = cumtrapz(temperature, specHeatMag./temperature);
entropyErr = nan(size(specHeatMag));
entropyErr(1) = (temperature(2)-temperature(1))/2/temperature(1)*specHeatErr(1); % See 04/11/2023 526 notes
for i = 2:length(specHeatMag)-1
    arg = nan(i, 1);
    arg(1) = ((temperature(2)-temperature(1))/2/temperature(1)*specHeatErr(1)).^2;
    for j = 2:i
        arg(j) = ((temperature(j+1)-temperature(j-1))/2/temperature(j)).^2.*specHeatErr(j).^2;
    end
    entropyErr(i) = sqrt(sum(arg));
end
entropyErr(end) = sqrt(sum([arg; ((temperature(end)-temperature(end-1))/2/temperature(end).*specHeatErr(end)).^2])); % See 04/11/2023 526 notes

disp(['Debye Temperature (K): ', num2str(x(1))])
disp(['Debye Temperature Uncertainty (K): ', num2str(xErr(1))])
disp(['Alpha: ', num2str(x(2))])
disp(['Alpha Uncertainty: ', num2str(xErr(2))])
disp(['Reduced Chi-Squared: ', num2str(fval)])
disp(['Threshold Reduced Chi-Squared: ', num2str(chiUpper)])

% Plot the observed and calculated specific heat
ha = zeros(3,1);
figure('Units', 'inches', 'Position', [0, 1.0, 3.375, 4.96])
t = tiledlayout(3,1);
temperaturePlot = linspace(min(temperature), max(temperature), 5e2);
specHeatPlot = deb(temperaturePlot, x(1));
ha(1) = nexttile;
hold on
ylabel('\itC\rm/\itT\rm (J/mol\cdotK^2)')
errorbar(temperature, specHeat./temperature, specHeatErr./temperature, 'o', 'MarkerFaceColor', 'w')
p1 = plot(temperaturePlot, specHeatPlot./temperaturePlot, 'r', 'LineWidth', 2);
legend(p1, 'Model', 'fontsize', 10, 'AutoUpdate', 'off')
set(gca, 'TickLength', [0.02, 0.025])
pbaspect([16 9 1])
box on
hold off

% Plot the magnetic contribution to the specific heat
ha(2) = nexttile;
hold on
ylabel('\Delta\itC\rm/\itT\rm (J/mol\cdotK^2 )')
errorbar(temperature, specHeatMag./temperature, specHeatErr./temperature, 'o', 'MarkerFaceColor', 'w')
set(gca, 'TickLength', [0.02, 0.025])
pbaspect([16 9 1])
box on
hold off

% Plot the magnetic entropy
ha(3) = nexttile;
hold on
xlabel('\itT\rm (K)')
ylabel('\itS\rm_{mag} (J/mol\cdotK)')
errorbar(temperature, entropy, entropyErr, 'o', 'MarkerFaceColor', 'w')
set(gca, 'TickLength', [0.02, 0.025])
pbaspect([16 9 1])
box on
hold off
t.TileSpacing = 'none';
t.Padding = 'compact';

% Plot uncertainty
ha2 = zeros(2,1);
figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
t2 = tiledlayout(2,1);
ha2(1) = nexttile;
hold on
title(['\Theta_D: ', num2str(round(x(1), 3)), '\pm', num2str(round(xErr(1), 3))])
xlabel('\Theta_D (K)')
ylabel('\chi^2_{r}')
scatter(paramTrial(:, 1), chiTrial(:, 1))
plot([paramTrial(interpPts(1, 1), 1), paramTrial(interpPts(2, 1), 1)], [chiTrial(interpPts(1, 1), 1), chiTrial(interpPts(2, 1), 1)], 'b')
plot([paramTrial(interpPts(3, 1), 1), paramTrial(interpPts(4, 1), 1)], [chiTrial(interpPts(3, 1), 1), chiTrial(interpPts(4, 1), 1)], 'b')
yline(chiUpper, 'Color', 'r', 'LineWidth', 3.0)
plot([paramLower(1), paramUpper(1)], [chiUpper, chiUpper], 'k-.o', 'LineWidth', 2.0)
xlim([-inf inf])
ylim([-inf chiUpper+4*(chiUpper-fval)])
set(gca,'FontSize',12)
box on
hold off

ha2(2) = nexttile;
hold on
title(['Alpha: ', num2str(round(x(2), 6)), '\pm', num2str(round(xErr(2), 6))])
xlabel('Alpha')
ylabel('\chi^2_{r}')
scatter(paramTrial(:, 2), chiTrial(:, 2))
plot([paramTrial(interpPts(1, 2), 2), paramTrial(interpPts(2, 2), 2)], [chiTrial(interpPts(1, 2), 2), chiTrial(interpPts(2, 2), 2)], 'b')
plot([paramTrial(interpPts(3, 2), 2), paramTrial(interpPts(4, 2), 2)], [chiTrial(interpPts(3, 2), 2), chiTrial(interpPts(4, 2), 2)], 'b')
yline(chiUpper, 'Color', 'r', 'LineWidth', 3.0)
plot([paramLower(2), paramUpper(2)], [chiUpper, chiUpper], 'k-.o', 'LineWidth', 2.0)
xlim([-inf inf])
ylim([fval-0.1 chiUpper+(chiUpper-fval)])
set(gca,'FontSize',12)
box on
hold off

% Plot the overall hc
figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.5])
hold on
ylabel('\itC\rm (J/mol\cdotK)')
errorbar(temperature, specHeat, specHeatErr, 'o', 'MarkerFaceColor', 'w')
p1 = plot(temperaturePlot, specHeatPlot, 'r', 'LineWidth', 2);
legend(p1, 'Model', 'fontsize', 10, 'AutoUpdate', 'off', 'Location', 'northwest')
set(gca, 'TickLength', [0.02, 0.025])
pbaspect([16 9 1])
box on
hold off

end