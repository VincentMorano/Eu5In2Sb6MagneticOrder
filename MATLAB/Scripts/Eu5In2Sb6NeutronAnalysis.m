%% Vincent Morano

% BT7 experiment on Eu5In2Sb6: 12/11/2020

% Notes on filenames:
% 854081:854135 are 1.6K, (H0L)
% 854138:854187 are 1.6K, (0KL)
% 854188:854203 are 18K, (0KL)
% 854204:854215 are 18K, (H0L)
% 854315:854323 are 18K, (H0L) positive
% 854324:854334 are 18K, (H0L) other quadrants
% 854335:854342 are 18K, (0KL) positive
% 854343:854350 are 18K, (0KL) other quadrants

% Exclude 4197 and 4339 since tried a scan in the wrong plane.

%% Fit strong nuclear 0KL sample peaks for FWHM versus Q

clear
close all

factor = 1e3; % Scale the observed intensities by this arbitrary prefactor
headerLines = 44;
bgThresh = inf; % Don't save points fit with bg higher than this.
alumLines = [38.1675, 44.3615, 64.5397, 77.5203]; % A4 values of alum lines
dropRange = [2.5, 2.5, 2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0, 0.2]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.1, 0.1]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\';
fileNum = [854189, 854191, 854194, 854196, 854203];
fileNames = strcat('fpx',string(fileNum),'.bt7');
scans = importDataBT7(fileLoc,fileNames,headerLines);

% Import data
for i = 1:length(scans)
    scans(i).fileNum = fileNum(i);
    
    if ismember(scans(i).fileNum, [854188:854203, 854335:854350])
        scans(i).type = 1;
    elseif ismember(scans(i).fileNum, 854138:854187)
        scans(i).type = 2;
    end
end

% Fit for integrated intensities
numParam = 4; % Constraining std. dev. to mean of prominent peaks.
model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./x(4)).*exp(-4.*log(2).*(a3-x(3)).^2./x(4).^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
for i = 1:length(scans)
    % Apply a scale factor to the intensities if desired
    scans(i).intMon = factor*scans(i).intMon;
    scans(i).intMonErr = factor*scans(i).intMonErr;
    scans(i).intTime = factor*scans(i).intTime;
    scans(i).intTimeErr = factor*scans(i).intTimeErr;
    
    scans(i).numParam = numParam;
    scans(i).includePowder = sum(abs(mean(scans(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scans(i).a3);
    
    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scans(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intMon-bg0))./sum(scans(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    fwhm0 = 0.3; % Assume close to 0.3 deg.
    area0 = (max(scans(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0; area0; cent0; fwhm0]; % Starting point for fitting.
    bglb = min(scans(i).intMon);
    bgub = max(scans(i).intMon);
    centerlb = scans(i).a3(3);
    centerub = scans(i).a3(end-2);
    fwhmlb = 0.2;
    fwhmub = 0.6;
    arealb = -(max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    areaub = (max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    lb = [bglb, arealb, centerlb, fwhmlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub, fwhmub]; % fmincon upper bounds on fitted variables
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2ErrCon(scans(i).intMon,scans(i).intMonErr,modelInput,scans(i).x0,errPts,fact,offset,lb,ub);
    scans(i).include=(scans(i).xFit(1)<bgThresh) & scans(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 200)';
    intCal = model(scans(i).xFit, a3Cal);
    
    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scans(i).meanHKL), '\right)$, ', num2str(scans(i).meanT, 5), ' K, ', num2str(scans(i).meanE, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(scans)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scans(i).a3, scans(i).intMon, scans(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'NucFWHMCurve.png'])
    
    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    for j = 1:length(scans(i).x0)
        subplot(2, 2, j)
        hold on
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(2, 2, 1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '$\pm$', num2str(round(scans(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '$\pm$', num2str(round(scans(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '$\pm$', num2str(round(scans(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 4)
    hold on
    title(['Gaussian FWHM: ', num2str(round(scans(i).xFit(4), 4)), '$\pm$', num2str(round(scans(i).xErr(4), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian FWHM (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'NucFWHMError.png'])
    pause(0.1)
end

% Calculate structure factors
for i = 1:length(scans)
    scans(i).W = 0;
    scans(i).Exp=Eu5In2Sb6121120;
    scans(i).lattice = GetLattice(scans(i).Exp);
    scans(i).a = scans(i).lattice.a;
    scans(i).b = scans(i).lattice.b;
    scans(i).c = scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).meanH.^2/scans(i).a.^2+scans(i).meanK.^2/scans(i).b.^2+scans(i).meanL.^2/scans(i).c.^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
    [scans(i).R0,scans(i).M] = ResMat(scans(i).Q,scans(i).W,scans(i).Exp);
    [scans(i).f2,scans(i).f2Err] = structFact(scans(i).R0, scans(i).M, scans(i).xFit(2), scans(i).xErr(2), scans(i).Q);
end

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('Nuclear FWHM Versus Q', 'Interpreter', 'LaTeX')
ylabel('$\Delta \mathrm{A}_3$', 'Interpreter', 'LaTeX')
xlabel('Q $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scans), arrayfun(@(x) x.xFit(4), scans), arrayfun(@(x) x.xErr(4), scans), 'o')
yline(mean(arrayfun(@(x) x.xFit(4), scans)), '-', ['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scans)))])
box on
hold off
saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\NucFWHM.png')

disp(['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scans)))])

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\NucFWHMStrFct.mat') % Save the Workspace

%% Fit strong k=0,0,1/2 magnetic 0KL sample peaks for FWHM versus Q

clear
close all

factor = 1e3; % Scale the observed intensities by this arbitrary prefactor
headerLines = 44;
bgThresh = inf; % Don't save points fit with bg higher than this.
alumLines = [38.1675, 44.3615, 64.5397, 77.5203]; % A4 values of alum lines
dropRange = [2.5, 2.5, 2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0, 0.2]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.8, 0.1]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\';
fileNum = [854139, 854145, 854147, 854150, 854164, 854171];
fileNames = strcat('fpx',string(fileNum),'.bt7');
scans = importDataBT7(fileLoc,fileNames,headerLines);

% Import data
for i = 1:length(scans)
    scans(i).fileNum = fileNum(i);
    
    if ismember(scans(i).fileNum, [854188:854203, 854335:854350])
        scans(i).type = 1;
    elseif ismember(scans(i).fileNum, 854138:854187)
        scans(i).type = 2;
    end
end

% Fit for integrated intensities
numParam = 4; % Constraining std. dev. to mean of prominent peaks.
model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./x(4)).*exp(-4.*log(2).*(a3-x(3)).^2./x(4).^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
for i = 1:length(scans)
    % Apply a scale factor to the intensities if desired
    scans(i).intMon = factor*scans(i).intMon;
    scans(i).intMonErr = factor*scans(i).intMonErr;
    scans(i).intTime = factor*scans(i).intTime;
    scans(i).intTimeErr = factor*scans(i).intTimeErr;
    
    scans(i).numParam = numParam;
    scans(i).includePowder = sum(abs(mean(scans(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scans(i).a3);
    
    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scans(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intMon-bg0))./sum(scans(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    fwhm0 = 0.3; % Assume close to 0.3 deg.
    area0 = (max(scans(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0; area0; cent0; fwhm0]; % Starting point for fitting.
    bglb = min(scans(i).intMon);
    bgub = max(scans(i).intMon);
    centerlb = scans(i).a3(3);
    centerub = scans(i).a3(end-2);
    fwhmlb = 0.2;
    fwhmub = 0.6;
    arealb = -(max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    areaub = (max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    lb = [bglb, arealb, centerlb, fwhmlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub, fwhmub]; % fmincon upper bounds on fitted variables
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2ErrCon(scans(i).intMon,scans(i).intMonErr,modelInput,scans(i).x0,errPts,fact,offset,lb,ub);
    scans(i).include=(scans(i).xFit(1)<bgThresh) & scans(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 200)';
    intCal = model(scans(i).xFit, a3Cal);
    
    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scans(i).meanHKL), '\right)$, ', num2str(scans(i).meanT, 5), ' K, ', num2str(scans(i).meanE, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(scans)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scans(i).a3, scans(i).intMon, scans(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'MagFWHMCurve.png'])
    
    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    for j = 1:length(scans(i).x0)
        subplot(2, 2, j)
        hold on
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(2, 2, 1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '$\pm$', num2str(round(scans(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '$\pm$', num2str(round(scans(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '$\pm$', num2str(round(scans(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 4)
    hold on
    title(['Gaussian FWHM: ', num2str(round(scans(i).xFit(4), 4)), '$\pm$', num2str(round(scans(i).xErr(4), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian FWHM (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'MagFWHMError.png'])
    pause(0.1)
end

% Calculate structure factors
for i = 1:length(scans)
    scans(i).W = 0;
    scans(i).Exp=Eu5In2Sb6121120;
    scans(i).lattice = GetLattice(scans(i).Exp);
    scans(i).a = scans(i).lattice.a;
    scans(i).b = scans(i).lattice.b;
    scans(i).c = scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).meanH.^2/scans(i).a.^2+scans(i).meanK.^2/scans(i).b.^2+scans(i).meanL.^2/scans(i).c.^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
    [scans(i).R0,scans(i).M] = ResMat(scans(i).Q,scans(i).W,scans(i).Exp);
    [scans(i).f2,scans(i).f2Err] = structFact(scans(i).R0, scans(i).M, scans(i).xFit(2), scans(i).xErr(2), scans(i).Q);
end

% Fit to the function for my magnetic fwhm
modelFWHM = @(x, Q) sqrt(x(1).^2 + (x(2)./Q).^2);
x0 = [0.3; 0.1];
errPts = 1e3;
fact = [0.2, 0.2];
offset = [0, 0];
[xFit,redChi2Fit,xErr,chiUpper,chiTrial,paramTrial,interpPts,slopes,intercepts,paramLower,paramUpper] = fitRedChi2Err(arrayfun(@(x) x.xFit(4), scans),arrayfun(@(x) x.xErr(4), scans),@(x) modelFWHM(x, arrayfun(@(x) x.Q, scans)),x0,errPts,fact,offset);

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('k=(0,0,1/2) FWHM Versus Q', 'Interpreter', 'LaTeX')
ylabel('$\Delta \mathrm{A}_3$', 'Interpreter', 'LaTeX')
xlabel('Q $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scans), arrayfun(@(x) x.xFit(4), scans), arrayfun(@(x) x.xErr(4), scans), 'o')
QPlot = linspace(min(arrayfun(@(x) x.Q, scans)), max(arrayfun(@(x) x.Q, scans)), 5e2);
plot(QPlot, modelFWHM(xFit, QPlot))
box on
hold off
saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\MagFWHM.png')

disp(['Delta A3r: ', num2str(xFit(1))])
disp(['Delta Q: ', num2str(xFit(2))])

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\MagFWHMStrFct.mat') % Save the Workspace

%% Fit for integrated intensities of 0KL sample rocking scans, calculate structure factors.
clear
close all

fwhmNuc = 0.38748;
fwhmMag = 0.39789;
factor = 1e3; % Scale the observed intensities by this arbitrary prefactor
headerLines = 44;
bgThresh = inf; % Don't save points fit with bg higher than this.
alumLines = [38.1675, 44.3615, 64.5397, 77.5203]; % A4 values of alum lines
dropRange = [2.5, 2.5, 2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.4]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\';
fileNum = [854138:854196, 854198:854203];
fileNames = strcat('fpx',string(fileNum),'.bt7');
scans = importDataBT7(fileLoc,fileNames,headerLines);

% Import data
for i = 1:length(scans)
    scans(i).fileNum = fileNum(i);
    
    if ismember(scans(i).fileNum, 854188:854203)
        scans(i).type = 1;
    elseif ismember(scans(i).fileNum, 854138:854187)
        scans(i).type = 2;
    end
end

% Fit for integrated intensities
numParam = 3; % Constraining std. dev. to mean of prominent peaks.
for i = 1:length(scans)
    if scans(i).type==1
        model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./fwhmNuc).*exp(-4.*log(2).*(a3-x(3)).^2./fwhmNuc.^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
        fwhm0 = fwhmNuc;
    elseif scans(i).type==2
        model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./fwhmMag).*exp(-4.*log(2).*(a3-x(3)).^2./fwhmMag.^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
        fwhm0 = fwhmMag;
    end
    
    % Apply a scale factor to the intensities if desired
    scans(i).intMon = factor*scans(i).intMon;
    scans(i).intMonErr = factor*scans(i).intMonErr;
    scans(i).intTime = factor*scans(i).intTime;
    scans(i).intTimeErr = factor*scans(i).intTimeErr;
    
    scans(i).numParam = numParam;
    scans(i).includePowder = sum(abs(mean(scans(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scans(i).a3);
    
    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scans(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intMon-bg0))./sum(scans(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    area0 = (max(scans(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0; area0; cent0]; % Starting point for fitting.
    bglb = min(scans(i).intMon);
    bgub = max(scans(i).intMon);
    centerlb = scans(i).a3(7);
    centerub = scans(i).a3(end-6);
    arealb = -(max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*1.4*fwhm0/2; % Raising fwhm by a factor
    areaub = (max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*1.4*fwhm0/2;
    lb = [bglb, arealb, centerlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub]; % fmincon upper bounds on fitted variables
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2ErrCon(scans(i).intMon,scans(i).intMonErr,modelInput,scans(i).x0,errPts,fact,offset,lb,ub);
    scans(i).include=(scans(i).xFit(1)<bgThresh) & scans(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 200)';
    intCal = model(scans(i).xFit, a3Cal);
    
    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scans(i).meanHKL), '\right)$, ', num2str(scans(i).meanT, 5), ' K, ', num2str(scans(i).meanE, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(scans)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scans(i).a3, scans(i).intMon, scans(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'Curve.png'])
    
    figure('Units', 'normalized', 'Position', [0.5, 0.1, 0.5, 0.8])
    for j = 1:length(scans(i).x0)
        subplot(3, 1, j)
        hold on
        set(gca, 'PlotBoxAspectRatio',[2 1 1]);
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(3, 1, 1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '$\pm$', num2str(round(scans(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(3, 1, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '$\pm$', num2str(round(scans(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(3, 1, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '$\pm$', num2str(round(scans(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'Error.png'])
    pause(0.1)
end

% Manually set to use reflections near powder lines if they look okay
scans(59).include = true;
scans(60).include = true;

% Calculate structure factors
for i = 1:length(scans)
    scans(i).W = 0;
    scans(i).Exp=Eu5In2Sb6121120;
    scans(i).lattice = GetLattice(scans(i).Exp);
    scans(i).a = scans(i).lattice.a;
    scans(i).b = scans(i).lattice.b;
    scans(i).c = scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).meanH.^2/scans(i).a.^2+scans(i).meanK.^2/scans(i).b.^2+scans(i).meanL.^2/scans(i).c.^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
    [scans(i).R0,scans(i).M] = ResMat(scans(i).Q,scans(i).W,scans(i).Exp);
    [scans(i).f2,scans(i).f2Err] = structFact(scans(i).R0, scans(i).M, scans(i).xFit(2), scans(i).xErr(2), scans(i).Q);
end
save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\BT7StrFct.mat') % Save the Workspace

%% Save relevant numbers to output file for Mantid to generate absorption-corrected .int files.
% Save to Excel workbook. Sort based on measurement set.
expName=["18K0KL" "1p6K0KL"];
for i=1:length(expName)
    ind=arrayfun(@(x) x.type, scans)==i.*ones(size(scans)) & arrayfun(@(x) x.include, scans); % Pick out scans of a certain set and only if bg is reasonable, avoid Al lines
    H=arrayfun(@(x) x.meanH, scans(ind)); % Convert to the desired arrays since they're cells.
    K=arrayfun(@(x) x.meanK, scans(ind));
    L=arrayfun(@(x) x.meanL, scans(ind));
    a4=arrayfun(@(x) x.meanA4, scans(ind));
    a3=arrayfun(@(x) x.xFit(3), scans(ind));
    f2=arrayfun(@(x) x.f2, scans(ind));
    f2Err=arrayfun(@(x) x.f2Err, scans(ind));
    fileNamesExcel=strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\', expName(i),'.xlsx');
    delete(fileNamesExcel)
    T=table(H', K', L', a4', a3', f2', f2Err', 'VariableNames',{'H', 'K', 'L', 'TwoTheta', 'Phi', 'F2', 'F2Error'}); % Using the fitted center for A3, averaging a4 values for 2-theta (should be identical anyway).
    writetable(T, fileNamesExcel);
end
%% Save structure factors to .int files.
load('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\BT7StrFct.mat') % Open the Workspace

fileID18K0KL=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb618K0KL.int','w');
formatSpec='Eu5In2Sb6 BT7 INT File: 18K, 35meV, (0KL)\n';
fprintf(fileID18K0KL, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID18K0KL, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID18K0KL, formatSpec);
fileID1p6K0KL000p5=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb61p6K0KL000p5.int','w');
formatSpec='Eu5In2Sb6 BT7 INT File: 1.6K, 35meV, (0KL) k=(0,0,1/2) Peaks\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='(4i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='1\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='1 0.0 0.0 0.5\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
for i=1:length(scans)
    if scans(i).type==1
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).meanH==floor(scans(i).meanH) && scans(i).meanK==floor(scans(i).meanK) && scans(i).meanL==floor(scans(i).meanL) && scans(i).include
            fprintf(fileID18K0KL, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==2
        formatSpec='%5i%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).meanL~=floor(scans(i).meanL) && scans(i).include
            fprintf(fileID1p6K0KL000p5, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL-0.5, 1, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    end
end
fclose(fileID18K0KL);
fclose(fileID1p6K0KL000p5);

%% Save spherical absorption corrected structure factors to .int files. Only valid for (0KL) sample.
load('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\BT7StrFct.mat') % Open the Workspace

% Determine the absorption factors
volSamp = 2; % Volume of sample in mm3
atomsEu = 10; % Per unit cell
atomsIn = 4;
vol = 843.684293; % Of unit cell in angstroms
energy = 35; % Energy in meV

sigAbsEu = 4530*sqrt(25.3/energy); % barns. Times 10^-22 for mm^2
sigAbsIn = 193.8*sqrt(25.3/energy);
densEu = atomsEu/vol; % Inverse cubic Ang. Times 10^21 for mm^-3
densIn = atomsIn/vol;
muAbs = (densEu*sigAbsEu + densIn*sigAbsIn)/10; % Absorption coefficient in inverse mm. 
radEff = (volSamp*3/4/pi)^(1/3); % Effective radius for spherical absorption correction in mm

fun=@(R, theta, mu, r, alpha, phi) sin(alpha).*(3./(4.*pi.*R.^3)).*exp(-mu.*(sqrt(R.^2-r.^2.*cos(alpha).^2-r.^2.*sin(alpha).^2.*sin(theta.*(pi/180)+phi).^2)+sqrt(R.^2-r.^2.*cos(alpha).^2-r.^2.*sin(alpha).^2.*sin(theta.*(pi/180)-phi).^2)-2.*r.*sin(theta.*(pi/180)).*sin(alpha).*sin(phi))).*r.^2; % Picked up the sin(alpha) because want to integrate over alpha rather than cos(alpha), then flipped integral bounds.
funA=@(R, theta, mu) integral3(@(r, alpha, phi) fun(R, theta, mu, r, alpha, phi), 0, R, 0, pi, 0, 2*pi); % Transmission coefficient. integral3 seems to plug in an array of values so need dots in front of operations in function. Mentioned in documentation.

AStar=zeros(length(scans), 1);
for i=1:length(scans)
    AStar(i)=1./funA(radEff, scans(i).meanA4/2, muAbs);
end

% Save the files
fileID18K0KL=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb618K0KLCorr.int','w');
formatSpec='Eu5In2Sb6 BT7 INT File: 18K, 35meV, (0KL), Corrected\n';
fprintf(fileID18K0KL, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID18K0KL, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID18K0KL, formatSpec);
fileID1p6K0KL000p5=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb61p6K0KL000p5Corr.int','w');
formatSpec='Eu5In2Sb6 BT7 INT File: 1.6K, 35meV, (0KL) k=(0,0,1/2) Peaks, Corrected\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='(4i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='1.5289 0 0\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='1\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
formatSpec='1 0.0 0.0 0.5\n';
fprintf(fileID1p6K0KL000p5, formatSpec);
for i=1:length(scans)
    if scans(i).type==1
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).meanH==floor(scans(i).meanH) && scans(i).meanK==floor(scans(i).meanK) && scans(i).meanL==floor(scans(i).meanL) && scans(i).include
            fprintf(fileID18K0KL, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==2
        formatSpec='%5i%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).meanL~=floor(scans(i).meanL) && scans(i).include
            fprintf(fileID1p6K0KL000p5, formatSpec, scans(i).meanH, scans(i).meanK, scans(i).meanL-0.5, 1, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    end
end
fclose(fileID18K0KL);
fclose(fileID1p6K0KL000p5);

%% Analysis for earlier SPINS experiment: 09/29/2020

% Sample is in (HK0)

% Notes on fileNamess:

% 1.5K | 03024:03029, 03031:03039, 03064, 03066, 03068:03070, 03072:03076, 
% 03087:03093

% 10K | 03040:03049, 03065, 03094:03098

% 20K | 03050:03052, 03055:03056, 03058:03063, 03077:03083, 03099:03105
%% Fit to the fwhm of the PM SPINS peaks
clear
close all

factor = 1e5; % Scale the observed intensities by this arbitrary prefactor
bgThresh=inf; % Don't save points fit with bg higher than this.
alumLines=[119.768, 174.449]; % A4 values of alum lines for 35meV
dropRange=[2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0, 0.2]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.8, 0.1]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\';
fileNum = [03050, 03052, 03055, 03060:03062, 03082];
fileNames = string(strcat('fpx',num2str(fileNum(:), '%05d'),'.ng5'))';
T = 20.*ones(size(fileNum));
mon = 700000.*ones(size(fileNum));
H = [0, 4, 3, 1, 4, 4, 2];
K = [4, 0, 3, 4, 1, 2, 3];
L = [0, 0, 0, 0, 0, 0, 0];
E = 5.*ones(size(fileNum));

% Import data
scans = importDataSPINS(fileLoc,fileNames,T,mon,H,K,L,E);
for i=1:length(fileNames)
    scans(i).fileNum = fileNum(i);
    scans(i).Exp=Eu5In2Sb6100820;
    scans(i).lattice=GetLattice(scans(i).Exp);
    scans(i).a=scans(i).lattice.a;
    scans(i).b=scans(i).lattice.b;
    scans(i).c=scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).H.^2/scans(i).a.^2+scans(i).K.^2/scans(i).b.^2+scans(i).L.^2/scans(i).c.^2);
    scans(i).a4=2.*asind(scans(i).Q./2./(2*pi./4.045));
    
    if ismember(scans(i).fileNum, [03024:03029, 03031:03039, 03064, 03066, 03068:03070, 03072:03076, 03087:03093])
        scans(i).type=1;
    elseif ismember(scans(i).fileNum, [03040:03049, 03065, 03094:03098])
        scans(i).type=2;
    elseif ismember(scans(i).fileNum, [03050:03052, 03055:03056, 03058:03063, 03077:03083, 03099:030105])
        scans(i).type=3;
    end
end

% Fit for integrated intensities
numParam = 4; % Constraining std. dev. to mean of prominent peaks.
model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./x(4)).*exp(-4.*log(2).*(a3-x(3)).^2./x(4).^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
for i = 1:length(scans)
    % Apply a scale factor to the intensities if desired
    scans(i).intMon = factor*scans(i).intMon;
    scans(i).intMonErr = factor*scans(i).intMonErr;
    
    scans(i).numParam = numParam;
    scans(i).includePowder = sum(abs(mean(scans(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scans(i).a3);
    
    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scans(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intMon-bg0))./sum(scans(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    fwhm0 = 0.3; % Assume close to 0.3 deg.
    area0 = (max(scans(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0, area0, cent0, fwhm0]; % Starting point for fitting.
    bglb = min(scans(i).intMon);
    bgub = max(scans(i).intMon);
    centerlb = scans(i).a3(3);
    centerub = scans(i).a3(end-2);
    fwhmlb = 0.2;
    fwhmub = 0.6;
    arealb = -(max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    areaub = (max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    lb = [bglb, arealb, centerlb, fwhmlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub, fwhmub]; % fmincon upper bounds on fitted variables
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2ErrCon(scans(i).intMon,scans(i).intMonErr,modelInput,scans(i).x0,errPts,fact,offset,lb,ub);
    scans(i).include=(scans(i).xFit(1)<bgThresh) & scans(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 200)';
    intCal = model(scans(i).xFit, a3Cal);
    
    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scans(i).HKL), '\right)$, ', num2str(scans(i).T, 5), ' K, ', num2str(scans(i).E, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(scans)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scans(i).a3, scans(i).intMon, scans(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'NucFWHMCurveSPINS.png'])
    
    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    for j = 1:length(scans(i).x0)
        subplot(2, 2, j)
        hold on
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(2, 2, 1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '$\pm$', num2str(round(scans(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '$\pm$', num2str(round(scans(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '$\pm$', num2str(round(scans(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 4)
    hold on
    title(['Gaussian FWHM: ', num2str(round(scans(i).xFit(4), 4)), '$\pm$', num2str(round(scans(i).xErr(4), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian FWHM (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'NucFWHMErrorSPINS.png'])
    pause(0.1)
end

% Calculate structure factors
for i = 1:length(scans)
    scans(i).W = 0;
    [scans(i).R0,scans(i).M] = ResMat(scans(i).Q,scans(i).W,scans(i).Exp);
    [scans(i).f2,scans(i).f2Err] = structFact(scans(i).R0, scans(i).M, scans(i).xFit(2), scans(i).xErr(2), scans(i).Q);
end

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('Nuclear FWHM Versus Q', 'Interpreter', 'LaTeX')
ylabel('$\Delta \mathrm{A}_3$', 'Interpreter', 'LaTeX')
xlabel('Q $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scans), arrayfun(@(x) x.xFit(4), scans), arrayfun(@(x) x.xErr(4), scans), 'o')
yline(mean(arrayfun(@(x) x.xFit(4), scans)), '-', ['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scans)))])
box on
hold off
saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\NucFWHMSPINS.png')

disp(['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scans)))])

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\NucFWHMStrFctSPINS.mat') % Save the Workspace

%% Fit to the fwhm of the 10 K k=0 SPINS peaks
clear
close all

factor = 1e5; % Scale the observed intensities by this arbitrary prefactor
bgThresh=inf; % Don't save points fit with bg higher than this.
alumLines=[119.768, 174.449]; % A4 values of alum lines for 35meV
dropRange=[2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0, 0.2]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.8, 0.1]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\';
fileNum = [03043, 03055, 03046, 03059, 03047, 03060, 03065, 03063, 03094, 03079, 03097, 03082];
fileNames = string(strcat('fpx',num2str(fileNum(:), '%05d'),'.ng5'))';
T = [10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20];
H = [3, 3, 2, 2, 1, 1, 1, 1, 3, 3, 2, 2];
K = [3, 3, 1, 1, 4, 4, 0, 0, 1, 1, 3, 3];
mon = 700000.*ones(size(fileNum));
L = zeros(size(fileNum));
E = 5.*ones(size(fileNum));

% Import data
scans = importDataSPINS(fileLoc,fileNames,T,mon,H,K,L,E);
for i=1:length(fileNames)
    scans(i).fileNum = fileNum(i);
    scans(i).Exp=Eu5In2Sb6100820;
    scans(i).lattice=GetLattice(scans(i).Exp);
    scans(i).a=scans(i).lattice.a;
    scans(i).b=scans(i).lattice.b;
    scans(i).c=scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).H.^2/scans(i).a.^2+scans(i).K.^2/scans(i).b.^2+scans(i).L.^2/scans(i).c.^2);
    scans(i).a4=2.*asind(scans(i).Q./2./(2*pi./4.045));
    
    if ismember(scans(i).fileNum, [03024:03029, 03031:03039, 03064, 03066, 03068:03070, 03072:03076, 03087:03093])
        scans(i).type=1;
    elseif ismember(scans(i).fileNum, [03040:03049, 03065, 03094:03098])
        scans(i).type=2;
    elseif ismember(scans(i).fileNum, [03050:03052, 03055:03056, 03058:03063, 03077:03083, 03099:030105])
        scans(i).type=3;
    end
end

% Pick out different temperatures then subtract using existing code then
% run as usual. Append to existing cells.
temps=arrayfun(@(x) x.T, scans);
ind10KScans=temps==10;
ind20KScans=temps==20;

int10K=arrayfun(@(x) x.intMon, scans(ind10KScans), 'UniformOutput', false);
int10KErr=arrayfun(@(x) x.intMonErr, scans(ind10KScans), 'UniformOutput', false);
hkl10K=cell2mat(arrayfun(@(x) x.HKL, scans(ind10KScans), 'UniformOutput', false)');
int20K=arrayfun(@(x) x.intMon, scans(ind20KScans), 'UniformOutput', false);
int20KErr=arrayfun(@(x) x.intMonErr, scans(ind20KScans), 'UniformOutput', false);
hkl20K=cell2mat(arrayfun(@(x) x.HKL, scans(ind20KScans), 'UniformOutput', false)');

corrEle10K=find(ind10KScans);

[hkl10KSubt, i10KSubtMag, i10KSubtPM]=intersect(hkl10K, hkl20K, 'rows');

int10KSubt=mat2cell(cell2mat(int10K(i10KSubtMag))-cell2mat(int20K(i10KSubtPM)),length(int10K{1}),ones(length(i10KSubtPM),1));
int10KSubtErr=mat2cell(sqrt(cell2mat(int10KErr(i10KSubtMag)).^2+cell2mat(int20KErr(i10KSubtPM)).^2),length(int10K{1}),ones(length(i10KSubtPM),1));

for i=1:length(hkl10KSubt)
    scans(end+1) = scans(corrEle10K(i10KSubtMag(i)));
    scans(end).intMon=int10KSubt{i};
    scans(end).intMonErr=int10KSubtErr{i};
    scans(end).type=5;
    scans(end).corrEle=corrEle10K(i10KSubtMag(i));
end

% Fit for integrated intensities
numParam = 4; % Constraining std. dev. to mean of prominent peaks.
model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./x(4)).*exp(-4.*log(2).*(a3-x(3)).^2./x(4).^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
scansFit = scans(arrayfun(@(x) x.type, scans)==5);
for i = 1:length(scansFit)
    % Apply a scale factor to the intensities if desired
    scansFit(i).intMon = factor*scansFit(i).intMon;
    scansFit(i).intMonErr = factor*scansFit(i).intMonErr;

    scansFit(i).numParam = numParam;
    scansFit(i).includePowder = sum(abs(mean(scansFit(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scansFit(i).a3);

    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scansFit(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scansFit(i).a3.*(scansFit(i).intMon-bg0))./sum(scansFit(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    fwhm0 = 0.3; % Assume close to 0.3 deg.
    area0 = (max(scansFit(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scansFit(i).x0 = [bg0; area0; cent0; fwhm0]; % Starting point for fitting.
    bglb = min(scansFit(i).intMon);
    bgub = max(scansFit(i).intMon);
    centerlb = scansFit(i).a3(3);
    centerub = scansFit(i).a3(end-2);
    fwhmlb = 0.2;
    fwhmub = 0.6;
    arealb = -(max(scansFit(i).intMon+scansFit(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    areaub = (max(scansFit(i).intMon+scansFit(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    lb = [bglb, arealb, centerlb, fwhmlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub, fwhmub]; % fmincon upper bounds on fitted variables
    [scansFit(i).xFit,scansFit(i).redChi2Fit,scansFit(i).xErr,scansFit(i).chiUpper,scansFit(i).chiTrial,scansFit(i).paramTrial,scansFit(i).interpPts,scansFit(i).slopes,scansFit(i).intercepts,scansFit(i).paramLower,scansFit(i).paramUpper] = fitRedChi2ErrCon(scansFit(i).intMon,scansFit(i).intMonErr,modelInput,scansFit(i).x0,errPts,fact,offset,lb,ub);
    scansFit(i).include=(scansFit(i).xFit(1)<bgThresh) & scansFit(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scansFit(i).a3), max(scansFit(i).a3), 200)';
    intCal = model(scansFit(i).xFit, a3Cal);

    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scansFit(i).HKL), '\right)$, ', num2str(scansFit(i).T, 5), ' K, ', num2str(scansFit(i).E, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(hkl10KSubt)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scansFit(i).a3, scansFit(i).intMon, scansFit(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scansFit(i).a3), max(scansFit(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scansFit(i).fileNum), 'Mag10KFWHMCurveSPINS.png'])

    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    for j = 1:length(scansFit(i).x0)
        subplot(2, 2, j)
        hold on
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scansFit(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scansFit(i).paramTrial(:, j), scansFit(i).chiTrial(:, j))
            plot([scansFit(i).paramTrial(scansFit(i).interpPts(1, j), j), scansFit(i).paramTrial(scansFit(i).interpPts(2, j), j)], [scansFit(i).chiTrial(scansFit(i).interpPts(1, j), j), scansFit(i).chiTrial(scansFit(i).interpPts(2, j), j)], 'b')
            plot([scansFit(i).paramTrial(scansFit(i).interpPts(3, j), j), scansFit(i).paramTrial(scansFit(i).interpPts(4, j), j)], [scansFit(i).chiTrial(scansFit(i).interpPts(3, j), j), scansFit(i).chiTrial(scansFit(i).interpPts(4, j), j)], 'b')
            yline(scansFit(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scansFit(i).paramLower(j), scansFit(i).paramUpper(j)], [scansFit(i).chiUpper, scansFit(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scansFit(i).redChi2Fit-0.1 scansFit(i).chiUpper+(scansFit(i).chiUpper-scansFit(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(2, 2, 1)
    hold on
    title(['Background: ', num2str(round(scansFit(i).xFit(1), 4)), '$\pm$', num2str(round(scansFit(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scansFit(i).xFit(2), 4)), '$\pm$', num2str(round(scansFit(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scansFit(i).xFit(3), 4)), '$\pm$', num2str(round(scansFit(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 4)
    hold on
    title(['Gaussian FWHM: ', num2str(round(scansFit(i).xFit(4), 4)), '$\pm$', num2str(round(scansFit(i).xErr(4), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian FWHM (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scansFit(i).fileNum), 'Mag10KFWHMErrorSPINS.png'])
    pause(0.1)
end

% Fit to the function for my magnetic fwhm
modelFWHM = @(x, Q) sqrt(x(1).^2 + (x(2)./Q).^2);
x0 = [0.3; 0.1];
errPts = 1e3;
fact = [0.2, 0.2];
offset = [0, 0];
[xFit,redChi2Fit,xErr,chiUpper,chiTrial,paramTrial,interpPts,slopes,intercepts,paramLower,paramUpper] = fitRedChi2Err(arrayfun(@(x) x.xFit(4), scansFit(arrayfun(@(x) x.type, scansFit)==5)),arrayfun(@(x) x.xErr(4), scansFit(arrayfun(@(x) x.type, scansFit)==5)),@(x) modelFWHM(x, arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==5))),x0,errPts,fact,offset);

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('10 K k=(0,0,0) FWHM Versus Q', 'Interpreter', 'LaTeX')
ylabel('$\Delta \mathrm{A}_3$', 'Interpreter', 'LaTeX')
xlabel('Q $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==5)), arrayfun(@(x) x.xFit(4), scansFit(arrayfun(@(x) x.type, scansFit)==5)), arrayfun(@(x) x.xErr(4), scansFit(arrayfun(@(x) x.type, scansFit)==5)), 'o')
QPlot = linspace(min(arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==5))), max(arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==5))), 5e2);
plot(QPlot, modelFWHM(xFit, QPlot))
box on
hold off
saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\Mag10KFWHMSPINS.png')

disp(['Delta A3r: ', num2str(xFit(1))])
disp(['Delta Q: ', num2str(xFit(2))])

disp(['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scansFit)))])

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\Mag10KFWHMStrFctSPINS.mat') % Save the Workspace

%% Fit to the fwhm of the 1.5 K k=0 SPINS peaks
clear
close all

factor = 1e5; % Scale the observed intensities by this arbitrary prefactor
bgThresh=inf; % Don't save points fit with bg higher than this.
alumLines=[119.768, 174.449]; % A4 values of alum lines for 35meV
dropRange=[2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0, 0.2]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.8, 0.1]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\';
fileNum = [03033, 03055, 03036, 03059, 03037, 03060, 03064, 03063, 03072, 03079, 03075, 03082];
fileNames = string(strcat('fpx',num2str(fileNum(:), '%05d'),'.ng5'))';
T = [1.5, 20, 1.5, 20, 1.5, 20, 1.5, 20, 1.5, 20, 1.5, 20];
H = [3, 3, 2, 2, 1, 1, 1, 1, 3, 3, 2, 2];
K = [3, 3, 1, 1, 4, 4, 0, 0, 1, 1, 3, 3];
mon = 700000.*ones(size(fileNum));
L = zeros(size(fileNum));
E = 5.*ones(size(fileNum));

% Import data
scans = importDataSPINS(fileLoc,fileNames,T,mon,H,K,L,E);
for i=1:length(fileNames)
    scans(i).fileNum = fileNum(i);
    scans(i).Exp=Eu5In2Sb6100820;
    scans(i).lattice=GetLattice(scans(i).Exp);
    scans(i).a=scans(i).lattice.a;
    scans(i).b=scans(i).lattice.b;
    scans(i).c=scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).H.^2/scans(i).a.^2+scans(i).K.^2/scans(i).b.^2+scans(i).L.^2/scans(i).c.^2);
    scans(i).a4=2.*asind(scans(i).Q./2./(2*pi./4.045));
    
    if ismember(scans(i).fileNum, [03024:03029, 03031:03039, 03064, 03066, 03068:03070, 03072:03076, 03087:03093])
        scans(i).type=1;
    elseif ismember(scans(i).fileNum, [03040:03049, 03065, 03094:03098])
        scans(i).type=2;
    elseif ismember(scans(i).fileNum, [03050:03052, 03055:03056, 03058:03063, 03077:03083, 03099:030105])
        scans(i).type=3;
    end
end

% Pick out different temperatures then subtract using existing code then
% run as usual. Append to existing cells.
temps=arrayfun(@(x) x.T, scans);
ind1p5KScans=temps==1.5;
ind20KScans=temps==20;

int1p5K=arrayfun(@(x) x.intMon, scans(ind1p5KScans), 'UniformOutput', false);
int1p5KErr=arrayfun(@(x) x.intMonErr, scans(ind1p5KScans), 'UniformOutput', false);
hkl1p5K=cell2mat(arrayfun(@(x) x.HKL, scans(ind1p5KScans), 'UniformOutput', false)');
int20K=arrayfun(@(x) x.intMon, scans(ind20KScans), 'UniformOutput', false);
int20KErr=arrayfun(@(x) x.intMonErr, scans(ind20KScans), 'UniformOutput', false);
hkl20K=cell2mat(arrayfun(@(x) x.HKL, scans(ind20KScans), 'UniformOutput', false)');

corrEle1p5K=find(ind1p5KScans);

[hkl1p5KSubt, i1p5KSubtMag, i1p5KSubtPM]=intersect(hkl1p5K, hkl20K, 'rows');

int1p5KSubt=mat2cell(cell2mat(int1p5K(i1p5KSubtMag))-cell2mat(int20K(i1p5KSubtPM)),length(int1p5K{1}),ones(length(i1p5KSubtPM),1));
int1p5KSubtErr=mat2cell(sqrt(cell2mat(int1p5KErr(i1p5KSubtMag)).^2+cell2mat(int20KErr(i1p5KSubtPM)).^2),length(int1p5K{1}),ones(length(i1p5KSubtPM),1));

for i=1:length(hkl1p5KSubt)
    scans(end+1) = scans(corrEle1p5K(i1p5KSubtMag(i)));
    scans(end).intMon = int1p5KSubt{i};
    scans(end).intMonErr = int1p5KSubtErr{i};
    scans(end).type = 4;
    scans(end).corrEle = corrEle1p5K(i1p5KSubtMag(i));
end

% Fit for integrated intensities
numParam = 4; % Constraining std. dev. to mean of prominent peaks.
model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./x(4)).*exp(-4.*log(2).*(a3-x(3)).^2./x(4).^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
scansFit = scans(arrayfun(@(x) x.type, scans)==4);
for i = 1:length(scansFit)
    % Apply a scale factor to the intensities if desired
    scansFit(i).intMon = factor*scansFit(i).intMon;
    scansFit(i).intMonErr = factor*scansFit(i).intMonErr;

    scansFit(i).numParam = numParam;
    scansFit(i).includePowder = sum(abs(mean(scansFit(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scansFit(i).a3);

    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scansFit(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scansFit(i).a3.*(scansFit(i).intMon-bg0))./sum(scansFit(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    fwhm0 = 0.3; % Assume close to 0.3 deg.
    area0 = (max(scansFit(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scansFit(i).x0 = [bg0; area0; cent0; fwhm0]; % Starting point for fitting.
    bglb = min(scansFit(i).intMon);
    bgub = max(scansFit(i).intMon);
    centerlb = scansFit(i).a3(3);
    centerub = scansFit(i).a3(end-2);
    fwhmlb = 0.2;
    fwhmub = 0.6;
    arealb = -(max(scansFit(i).intMon+scansFit(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    areaub = (max(scansFit(i).intMon+scansFit(i).intMonErr)-bglb)*sqrt(pi/log(2))*fwhmub/2;
    lb = [bglb, arealb, centerlb, fwhmlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub, fwhmub]; % fmincon upper bounds on fitted variables
    [scansFit(i).xFit,scansFit(i).redChi2Fit,scansFit(i).xErr,scansFit(i).chiUpper,scansFit(i).chiTrial,scansFit(i).paramTrial,scansFit(i).interpPts,scansFit(i).slopes,scansFit(i).intercepts,scansFit(i).paramLower,scansFit(i).paramUpper] = fitRedChi2ErrCon(scansFit(i).intMon,scansFit(i).intMonErr,modelInput,scansFit(i).x0,errPts,fact,offset,lb,ub);
    scansFit(i).include=(scansFit(i).xFit(1)<bgThresh) & scansFit(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scansFit(i).a3), max(scansFit(i).a3), 200)';
    intCal = model(scansFit(i).xFit, a3Cal);

    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scansFit(i).HKL), '\right)$, ', num2str(scansFit(i).T, 5), ' K, ', num2str(scansFit(i).E, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(hkl1p5KSubt)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scansFit(i).a3, scansFit(i).intMon, scansFit(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scansFit(i).a3), max(scansFit(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scansFit(i).fileNum), 'Mag1p5KFWHMCurveSPINS.png'])

    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    for j = 1:length(scansFit(i).x0)
        subplot(2, 2, j)
        hold on
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scansFit(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scansFit(i).paramTrial(:, j), scansFit(i).chiTrial(:, j))
            plot([scansFit(i).paramTrial(scansFit(i).interpPts(1, j), j), scansFit(i).paramTrial(scansFit(i).interpPts(2, j), j)], [scansFit(i).chiTrial(scansFit(i).interpPts(1, j), j), scansFit(i).chiTrial(scansFit(i).interpPts(2, j), j)], 'b')
            plot([scansFit(i).paramTrial(scansFit(i).interpPts(3, j), j), scansFit(i).paramTrial(scansFit(i).interpPts(4, j), j)], [scansFit(i).chiTrial(scansFit(i).interpPts(3, j), j), scansFit(i).chiTrial(scansFit(i).interpPts(4, j), j)], 'b')
            yline(scansFit(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scansFit(i).paramLower(j), scansFit(i).paramUpper(j)], [scansFit(i).chiUpper, scansFit(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scansFit(i).redChi2Fit-0.1 scansFit(i).chiUpper+(scansFit(i).chiUpper-scansFit(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(2, 2, 1)
    hold on
    title(['Background: ', num2str(round(scansFit(i).xFit(1), 4)), '$\pm$', num2str(round(scansFit(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scansFit(i).xFit(2), 4)), '$\pm$', num2str(round(scansFit(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scansFit(i).xFit(3), 4)), '$\pm$', num2str(round(scansFit(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    subplot(2, 2, 4)
    hold on
    title(['Gaussian FWHM: ', num2str(round(scansFit(i).xFit(4), 4)), '$\pm$', num2str(round(scansFit(i).xErr(4), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian FWHM (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scansFit(i).fileNum), 'Mag1p5KFWHMErrorSPINS.png'])
    pause(0.1)
end

% Fit to the function for my magnetic fwhm
modelFWHM = @(x, Q) sqrt(x(1).^2 + (x(2)./Q).^2);
x0 = [0.3; 0.1];
errPts = 1e3;
fact = [0.2, 0.2];
offset = [0, 0];
[xFit,redChi2Fit,xErr,chiUpper,chiTrial,paramTrial,interpPts,slopes,intercepts,paramLower,paramUpper] = fitRedChi2Err(arrayfun(@(x) x.xFit(4), scansFit(arrayfun(@(x) x.type, scansFit)==4)),arrayfun(@(x) x.xErr(4), scansFit(arrayfun(@(x) x.type, scansFit)==4)),@(x) modelFWHM(x, arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==4))),x0,errPts,fact,offset);

figure('Units', 'normalized', 'Position', [0.25, 0.3, 0.5, 0.6])
hold on
title('1.5 K k=(0,0,0) FWHM Versus Q', 'Interpreter', 'LaTeX')
ylabel('$\Delta \mathrm{A}_3$', 'Interpreter', 'LaTeX')
xlabel('Q $(\mathrm{\AA}^{-1})$', 'Interpreter', 'LaTeX')
errorbar(arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==4)), arrayfun(@(x) x.xFit(4), scansFit(arrayfun(@(x) x.type, scansFit)==4)), arrayfun(@(x) x.xErr(4), scansFit(arrayfun(@(x) x.type, scansFit)==4)), 'o')
QPlot = linspace(min(arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==4))), max(arrayfun(@(x) x.Q, scansFit(arrayfun(@(x) x.type, scansFit)==4))), 5e2);
plot(QPlot, modelFWHM(xFit, QPlot))
box on
hold off
saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\Mag1p5KFWHMSPINS.png')

disp(['Delta A3r: ', num2str(xFit(1))])
disp(['Delta Q: ', num2str(xFit(2))])

disp(['Mean: ', num2str(mean(arrayfun(@(x) x.xFit(4), scansFit)))])

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\Mag1p5KFWHMStrFctSPINS.mat') % Save the Workspace

%% Fit and save SPINS files with the raw nuclear intensities subtracted off magnetic peaks
clear
close all

fwhmNuc = 0.48742;
fwhmMag10K = @(Q) sqrt(0.48387.^2 + (0.13958./Q).^2);
fwhmMag1p5K = @(Q) sqrt(0.4703.^2 + (0.18693./Q).^2);

factor = 1e4; % Scale the observed intensities by this arbitrary prefactor
headerLines = 44;
bgThresh=400; % Don't save points fit with bg higher than this.
alumLines=[119.768, 174.449]; % A4 values of alum lines for 35meV
dropRange=[2.5, 2.5]; % Angular range about alum lines to drop if within. +-
errPts = 5e2; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.5, 0.0, 0.2]; % Determines range of iterated values for errorbar. [BG, area, center, fwhm]
offset = [0.2, 0.5, 0.1, 0.1]; % Determines range of iterated values for errorbar

% Import settings. Set scan types in for loop below.
fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\';
fileNum = [03024:03029, 03031:03052, 03055:03056, 03058:03063, 03064, 03065, 03068:03070, 03072:03076, 03078:03083, 03087:03093, 03094:03098, 03099:03105];
fileNames = string(strcat('fpx',num2str(fileNum(:), '%05d'),'.ng5'))';
T = [1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,10,10,10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,20,20,20,1.5,10,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,20,20,20,20,20,20,1.5,1.5,1.5,1.5,1.5,1.5,1.5,10,10,10,10,10,20,20,20,20,20,20,20];
H = [0, 0, 0, 0, 2, 3, 3, 4, 3, 2, 1, 2, 1, 4, 4, 0, 3, 4, 3, 2, 1, 2, 1, 4, 4, 0, 3, 4, 3, 2, 1, 2, 1, 4, 4, 1, 1, 1, 2, 1/2, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2, 2, 3, 0, 0, 0, 1, 2, 2, 3, 3, 1, 2, 2, 3, 0, 0, 0, 1, 2, 2, 3];
K = [1, 2, 3, 4, 0, 0, 0, 0, 3, 2, 2, 1, 4, 1, 2, 4, 0, 0, 3, 2, 2, 1, 4, 1, 2, 4, 0, 0, 3, 2, 2, 1, 4, 1, 2, 0, 0, 0, 1/2, 0, 1, 1, 3, 4, 3, 2, 1, 1, 3, 4, 3, 2, 1, 2, 3, 1, 2, 0, 0, 1, 3, 4, 3, 2, 1, 2, 3, 1, 2, 0, 0];
mon = [repmat(356292, [1,6]), repmat(700000, [1,65])];
L = zeros(size(fileNum));
E = 5.*ones(size(fileNum));

% Import data
scans = importDataSPINS(fileLoc,fileNames,T,mon,H,K,L,E);
for i=1:length(fileNames)
    scans(i).fileNum = fileNum(i);
    scans(i).Exp=Eu5In2Sb6100820;
    scans(i).lattice=GetLattice(scans(i).Exp);
    scans(i).a=scans(i).lattice.a;
    scans(i).b=scans(i).lattice.b;
    scans(i).c=scans(i).lattice.c;
    scans(i).Q=2*pi*sqrt(scans(i).H.^2/scans(i).a.^2+scans(i).K.^2/scans(i).b.^2+scans(i).L.^2/scans(i).c.^2);
    scans(i).a4=2.*asind(scans(i).Q./2./(2*pi./4.045));
    
    if ismember(scans(i).fileNum, [03024:03029, 03031:03039, 03064, 03066, 03068:03070, 03072:03076, 03087:03093])
        scans(i).type=1;
    elseif ismember(scans(i).fileNum, [03040:03049, 03065, 03094:03098])
        scans(i).type=2;
    elseif ismember(scans(i).fileNum, [03050:03052, 03055:03056, 03058:03063, 03077:03083, 03099:030105])
        scans(i).type=3;
    end
end

% Pick out different temperatures then subtract using existing code then
% run as usual. Append to existing cells.
temps=arrayfun(@(x) x.T, scans);
ind1p5KScans=temps==1.5;
ind10KScans=temps==10;
ind20KScans=temps==20;

int1p5K=arrayfun(@(x) x.intMon, scans(ind1p5KScans), 'UniformOutput', false);
int1p5KErr=arrayfun(@(x) x.intMonErr, scans(ind1p5KScans), 'UniformOutput', false);
hkl1p5K=cell2mat(arrayfun(@(x) x.HKL, scans(ind1p5KScans), 'UniformOutput', false)');
int10K=arrayfun(@(x) x.intMon, scans(ind10KScans), 'UniformOutput', false);
int10KErr=arrayfun(@(x) x.intMonErr, scans(ind10KScans), 'UniformOutput', false);
hkl10K=cell2mat(arrayfun(@(x) x.HKL, scans(ind10KScans), 'UniformOutput', false)');
int20K=arrayfun(@(x) x.intMon, scans(ind20KScans), 'UniformOutput', false);
int20KErr=arrayfun(@(x) x.intMonErr, scans(ind20KScans), 'UniformOutput', false);
hkl20K=cell2mat(arrayfun(@(x) x.HKL, scans(ind20KScans), 'UniformOutput', false)');

corrEle1p5K=find(ind1p5KScans); % I've selected out the 1.5 K scans and would like their original indices from the scans structure
corrEle10K=find(ind10KScans);

[hkl1p5KSubt, i1p5KSubtMag, i1p5KSubtPM]=intersect(hkl1p5K, hkl20K, 'rows');
[hkl10KSubt, i10KSubtMag, i10KSubtPM]=intersect(hkl10K, hkl20K, 'rows');

int1p5KSubt=mat2cell(cell2mat(int1p5K(i1p5KSubtMag))-cell2mat(int20K(i1p5KSubtPM)),length(int1p5K{1}),ones(length(i1p5KSubtPM),1));
int1p5KSubtErr=mat2cell(sqrt(cell2mat(int1p5KErr(i1p5KSubtMag)).^2+cell2mat(int20KErr(i1p5KSubtPM)).^2),length(int1p5K{1}),ones(length(i1p5KSubtPM),1));
int10KSubt=mat2cell(cell2mat(int10K(i10KSubtMag))-cell2mat(int20K(i10KSubtPM)),length(int10K{1}),ones(length(i10KSubtPM),1));
int10KSubtErr=mat2cell(sqrt(cell2mat(int10KErr(i10KSubtMag)).^2+cell2mat(int20KErr(i10KSubtPM)).^2),length(int10K{1}),ones(length(i10KSubtPM),1));

for i=1:length(hkl1p5KSubt)
    scans(end+1) = scans(corrEle1p5K(i1p5KSubtMag(i))); % I've intersected the 1.5 K and 20 K arrays but want the original indices from the scans structure
    scans(end).intMon = int1p5KSubt{i};
    scans(end).intMonErr = int1p5KSubtErr{i};
    scans(end).type = 4;
    scans(end).corrEle = corrEle1p5K(i1p5KSubtMag(i));
end
for i=1:length(hkl10KSubt)
    scans(end+1) = scans(corrEle10K(i10KSubtMag(i)));
    scans(end).intMon = int10KSubt{i};
    scans(end).intMonErr = int10KSubtErr{i};
    scans(end).type = 5;
    scans(end).corrEle = corrEle10K(i10KSubtMag(i));
end

% Fit for integrated intensities
numParam = 3; % Constraining std. dev. to mean of prominent peaks.
for i = 1:length(scans)
    if scans(i).type==3
        model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./fwhmNuc).*exp(-4.*log(2).*(a3-x(3)).^2./fwhmNuc.^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
        fwhm0 = fwhmNuc;
    elseif scans(i).type==4 || scans(i).type==1
        fwhm0 = fwhmMag10K(scans(i).Q);
        model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./fwhm0).*exp(-4.*log(2).*(a3-x(3)).^2./fwhm0.^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
    elseif scans(i).type==5 || scans(i).type==2
        fwhm0 = fwhmMag1p5K(scans(i).Q);
        model = @(x, a3) x(1)+sqrt(log(2)./pi).*(2.*x(2)./fwhm0).*exp(-4.*log(2).*(a3-x(3)).^2./fwhm0.^2); % Gaussian peak in terms of area and fwhm rather than max height and standard deviation.
    end
    % Apply a scale factor to the intensities if desired
    scans(i).intMon = factor*scans(i).intMon;
    scans(i).intMonErr = factor*scans(i).intMonErr;
    
    scans(i).numParam = numParam;
    scans(i).includePowder = sum(abs(mean(scans(i).a4).*ones(size(alumLines))-alumLines)>dropRange)==length(alumLines); % Record if reflection is close to a powder line (0 if so). These were a3 scans so all a4 values should be identical.
    modelInput = @(x) model(x, scans(i).a3);
    
    % Set initial values, fit my minimizing reduced chi-squared
    tmp = sort(scans(i).intMon); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(scans(i).a3.*(scans(i).intMon-bg0))./sum(scans(i).intMon-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    area0 = (max(scans(i).intMon)-bg0)*sqrt(pi/log(2))*fwhm0/2; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    scans(i).x0 = [bg0; area0; cent0]; % Starting point for fitting.
    bglb = min(scans(i).intMon);
    bgub = max(scans(i).intMon);
    centerlb = scans(i).a3(8);
    centerub = scans(i).a3(end-7);
    arealb = -(max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*1.5*fwhm0/2;
    areaub = (max(scans(i).intMon+scans(i).intMonErr)-bglb)*sqrt(pi/log(2))*1.5*fwhm0/2;
    lb = [bglb, arealb, centerlb]; % fmincon lower bounds on fitted variables. [BG, area, center, fwhm]
    ub = [bgub, areaub, centerub]; % fmincon upper bounds on fitted variables
    [scans(i).xFit,scans(i).redChi2Fit,scans(i).xErr,scans(i).chiUpper,scans(i).chiTrial,scans(i).paramTrial,scans(i).interpPts,scans(i).slopes,scans(i).intercepts,scans(i).paramLower,scans(i).paramUpper] = fitRedChi2ErrCon(scans(i).intMon,scans(i).intMonErr,modelInput,scans(i).x0,errPts,fact,offset,lb,ub);
    scans(i).include=(scans(i).xFit(1)<bgThresh) & scans(i).includePowder;

    % Plot fits
    a3Cal = linspace(min(scans(i).a3), max(scans(i).a3), 200)';
    intCal = model(scans(i).xFit, a3Cal);
    
    close all
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title(['Rocking Scan Fit: $\left(', num2str(scans(i).HKL), '\right)$, ', num2str(scans(i).T, 5), ' K, ', num2str(scans(i).E, 3), ' meV, ', '$\left( \frac{', num2str(i), '}{', num2str(length(scans)), '} \right)$'], 'Interpreter', 'LaTeX')
    xlabel('A3 (deg.)', 'Interpreter', 'LaTeX')
    ylabel('Intensity $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    errorbar(scans(i).a3, scans(i).intMon, scans(i).intMonErr, 'o')
    plot(a3Cal, intCal)
    set(gca,'FontSize',12)
    xlim([min(scans(i).a3), max(scans(i).a3)])
    box on
    axis square
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'CurveSPINS.png'])
    
    figure('Units', 'normalized', 'Position', [0.5, 0.1, 0.5, 0.8])
    for j = 1:length(scans(i).x0)
        subplot(3, 1, j)
        hold on
        ylabel('$\chi^2_{\mathrm{r}}$', 'Interpreter', 'LaTeX')
        if ~any(isnan(scans(i).interpPts(:, j))) % Check to see if an errorbar was determined before plotting
            scatter(scans(i).paramTrial(:, j), scans(i).chiTrial(:, j))
            plot([scans(i).paramTrial(scans(i).interpPts(1, j), j), scans(i).paramTrial(scans(i).interpPts(2, j), j)], [scans(i).chiTrial(scans(i).interpPts(1, j), j), scans(i).chiTrial(scans(i).interpPts(2, j), j)], 'b')
            plot([scans(i).paramTrial(scans(i).interpPts(3, j), j), scans(i).paramTrial(scans(i).interpPts(4, j), j)], [scans(i).chiTrial(scans(i).interpPts(3, j), j), scans(i).chiTrial(scans(i).interpPts(4, j), j)], 'b')
            yline(scans(i).chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([scans(i).paramLower(j), scans(i).paramUpper(j)], [scans(i).chiUpper, scans(i).chiUpper], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([scans(i).redChi2Fit-0.1 scans(i).chiUpper+(scans(i).chiUpper-scans(i).redChi2Fit)])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(3, 1, 1)
    hold on
    title(['Background: ', num2str(round(scans(i).xFit(1), 4)), '$\pm$', num2str(round(scans(i).xErr(1), 4))], 'Interpreter', 'LaTeX')
    xlabel('Background $\left( \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(3, 1, 2)
    hold on
    title(['Integrated Intensity: ', num2str(round(scans(i).xFit(2), 4)), '$\pm$', num2str(round(scans(i).xErr(2), 4))], 'Interpreter', 'LaTeX')
    xlabel('Integrated Intensity $\left( \mathrm{deg.} \cdot \frac{\mbox{det cts}}{\mbox{mon cts}} \right)$', 'Interpreter', 'LaTeX')
    hold off
    subplot(3, 1, 3)
    hold on
    title(['Gaussian Center: ', num2str(round(scans(i).xFit(3), 4)), '$\pm$', num2str(round(scans(i).xErr(3), 4))], 'Interpreter', 'LaTeX')
    xlabel('Gaussian Center (deg.)', 'Interpreter', 'LaTeX')
    hold off
    saveas(gcf,['C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\StructureFactors\fpx', num2str(scans(i).fileNum), 'ErrorSPINS.png'])
    pause(0.1)
end

% Calculate structure factors
for i = 1:length(scans)
    scans(i).W = 0;
    [scans(i).R0,scans(i).M] = ResMat(scans(i).Q,scans(i).W,scans(i).Exp);
    [scans(i).f2,scans(i).f2Err] = structFact(scans(i).R0, scans(i).M, scans(i).xFit(2), scans(i).xErr(2), scans(i).Q);
end
save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\SPINSStrFct.mat') % Save the Workspace

%% Save relevant numbers to output file for Mantid to generate absorption-corrected .int files.
% Save to Excel workbook. Sort based on measurement set.
expName=["1p5KHK0" "10KHK0" "20KHK0" "1p5KHK0Subt" "10KHK0Subt"];
for i=1:length(expName)
    ind=arrayfun(@(x) x.type, scans)==i.*ones(size(scans)) & arrayfun(@(x) x.include, scans); % Pick out scans of a certain set and only if bg is reasonable, avoid Al lines
    H=arrayfun(@(x) x.H, scans(ind)); % Convert to the desired arrays since they're cells.
    K=arrayfun(@(x) x.K, scans(ind));
    L=arrayfun(@(x) x.L, scans(ind));
    a4=arrayfun(@(x) x.a4, scans(ind));
    a3=arrayfun(@(x) x.meanA3, scans(ind));
    f2=arrayfun(@(x) x.f2, scans(ind));
    f2Err=arrayfun(@(x) x.f2Err, scans(ind));
    fileNamesExcel=strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\', expName(i),'.xlsx');
    delete(fileNamesExcel)
    T=table(H', K', L', a4', a3', f2', f2Err', 'VariableNames',{'H', 'K', 'L', 'TwoTheta', 'Phi', 'F2', 'F2Error'}); % Using the fitted center for A3, averaging a4 values for 2-theta (should be identical anyway).
    writetable(T, fileNamesExcel);
end
%% Save structure factors to .int files.
fileID1p5KHK0=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb61p5KHK0.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 1.5K, 5meV, (HK0)\n';
fprintf(fileID1p5KHK0, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p5KHK0, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID1p5KHK0, formatSpec);
fileID10KHK0=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb610KHK0.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 10K, 5meV, (HK0)\n';
fprintf(fileID10KHK0, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID10KHK0, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID10KHK0, formatSpec);
fileID20KHK0=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb620KHK0.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 20K, 5meV, (HK0)\n';
fprintf(fileID20KHK0, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID20KHK0, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID20KHK0, formatSpec);
fileID1p5KHK0Subt=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb61p5KHK0Subt.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 1.5K, 5meV, (HK0), Subtracted\n';
fprintf(fileID1p5KHK0Subt, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p5KHK0Subt, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID1p5KHK0Subt, formatSpec);
fileID10KHK0Subt=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb610KHK0Subt.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 10K, 5meV, (HK0), Subtracted\n';
fprintf(fileID10KHK0Subt, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID10KHK0Subt, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID10KHK0Subt, formatSpec);
for i=1:length(scans)
    if scans(i).type==1
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID1p5KHK0, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==2
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID10KHK0, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==3
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID20KHK0, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==4
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID1p5KHK0Subt, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==5
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID10KHK0Subt, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2, scans(i).f2Err, 1, 0.00, 0.00, 0.00);
        end
    end
end
fclose(fileID1p5KHK0);
fclose(fileID10KHK0);
fclose(fileID20KHK0);
fclose(fileID1p5KHK0Subt);
fclose(fileID10KHK0Subt);

%% Save spherical absorption corrected structure factors to .int files.
load('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\SPINSStrFct.mat') % Open the Workspace

% Determine the absorption factors
volSamp = 12; % Volume of sample in mm3
atomsEu = 10; % Per unit cell
atomsIn = 4;
vol = 843.684293; % Of unit cell in angstroms
energy = 5; % Energy in meV

sigAbsEu = 4530*sqrt(25.3/energy); % barns. Times 10^-22 for mm^2
sigAbsIn = 193.8*sqrt(25.3/energy);
densEu = atomsEu/vol; % Inverse cubic Ang. Times 10^21 for mm^-3
densIn = atomsIn/vol;
muAbs = (densEu*sigAbsEu + densIn*sigAbsIn)/10; % Absorption coefficient in inverse mm. 
radEff = (volSamp*3/4/pi)^(1/3); % Effective radius for spherical absorption correction in mm

fun=@(R, theta, mu, r, alpha, phi) sin(alpha).*(3./(4.*pi.*R.^3)).*exp(-mu.*(sqrt(R.^2-r.^2.*cos(alpha).^2-r.^2.*sin(alpha).^2.*sin(theta.*(pi/180)+phi).^2)+sqrt(R.^2-r.^2.*cos(alpha).^2-r.^2.*sin(alpha).^2.*sin(theta.*(pi/180)-phi).^2)-2.*r.*sin(theta.*(pi/180)).*sin(alpha).*sin(phi))).*r.^2; % Picked up the sin(alpha) because want to integrate over alpha rather than cos(alpha), then flipped integral bounds.
funA=@(R, theta, mu) integral3(@(r, alpha, phi) fun(R, theta, mu, r, alpha, phi), 0, R, 0, pi, 0, 2*pi, 'AbsTol', 1e-8, 'RelTol', 1e-4); % Transmission coefficient. integral3 seems to plug in an array of values so need dots in front of operations in function. Mentioned in documentation.

AStar=zeros(length(scans), 1);
for i=1:length(scans)
    AStar(i)=1./funA(radEff, scans(i).a4/2, muAbs);
end

% Save the files
fileID1p5KHK0=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb61p5KHK0Corr.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 1.5K, 5meV, (HK0), Corrected\n';
fprintf(fileID1p5KHK0, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p5KHK0, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID1p5KHK0, formatSpec);
fileID10KHK0=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb610KHK0Corr.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 10K, 5meV, (HK0), Corrected\n';
fprintf(fileID10KHK0, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID10KHK0, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID10KHK0, formatSpec);
fileID20KHK0=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb620KHK0Corr.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 20K, 5meV, (HK0), Corrected\n';
fprintf(fileID20KHK0, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID20KHK0, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID20KHK0, formatSpec);
fileID1p5KHK0Subt=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb61p5KHK0SubtCorr.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 1.5K, 5meV, (HK0), Subtracted, Corrected\n';
fprintf(fileID1p5KHK0Subt, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID1p5KHK0Subt, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID1p5KHK0Subt, formatSpec);
fileID10KHK0Subt=fopen('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\Data\Eu5In2Sb6\Eu5In2Sb610KHK0SubtCorr.int','w');
formatSpec='Eu5In2Sb6 SPINS INT File: 10K, 5meV, (HK0), Subtracted, Corrected\n';
fprintf(fileID10KHK0Subt, formatSpec);
formatSpec='(3i5,2f15.6,i5,3f7.2)\n';
fprintf(fileID10KHK0Subt, formatSpec);
formatSpec='4.0450 0 0\n';
fprintf(fileID10KHK0Subt, formatSpec);
for i=1:length(scans)
    if scans(i).type==1
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID1p5KHK0, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==2
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID10KHK0, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==3
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID20KHK0, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==4
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID1p5KHK0Subt, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    elseif scans(i).type==5
        formatSpec='%5i%5i%5i%15.6f%15.6f%5i%7.2f%7.2f%7.2f\n';
        if scans(i).H==floor(scans(i).H) && scans(i).K==floor(scans(i).K) && scans(i).L==floor(scans(i).L) && scans(i).include
            fprintf(fileID10KHK0Subt, formatSpec, scans(i).H, scans(i).K, scans(i).L, scans(i).f2*AStar(i), scans(i).f2Err*AStar(i), 1, 0.00, 0.00, 0.00);
        end
    end
end
fclose(fileID1p5KHK0);
fclose(fileID10KHK0);
fclose(fileID20KHK0);
fclose(fileID1p5KHK0Subt);
fclose(fileID10KHK0Subt);
