%% Plots for the Eu5In2Sb6 manuscript figures. Vincent Morano, 05/11/2022
%% Figure 1
% Debye Fit: No Field 12/25/20 Without Threshold on Calculated HC
clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 10)

% Initial values for fitting
transition7K = 7.2028; % Midpoint of sharp edge
transition14K = 14.0443; % Midpoint of sharp edge
fitLimHigh = 200; % Don't fit to temperatures beyond this upper limit
fitLimLow = 90; % Don't fit to temperatures beyond this lower limit
debTemp0 = 170; % K
anharm0 = 1e-4; % Initial fitting value for a rescaling to match Dulong-Petit. Can think of it as an effective mass factor I suppose.
errPts = 1e2;
fact = [0, 0]; % Debye temperature, T-linear anharmonic correction
offset = [5, 0.00001];

% Constants
mass = 5.19; % Sample mass in mg
massErr = 0.03; % IQM scale uncertainty in mg, see 04/06/2023 526 notes
molarMass = 5*(151.96)+2*(114.82)+6*(121.76);
n = 13; % Atoms in formula unit
N = 6.022e23; % Avogadro's number
kB = 1.380649e-23; % Boltzmann constant
R = N*kB; % Gas constant in J/mol-K

fileLoc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\HeatCapacity\Eu5In2Sb6\HC12252020\';
fileName = 'hc_sample_12232020_1p8_300K.dat';
file = readcell(strcat(fileLoc, fileName), 'FileType', 'text', 'NumHeaderLines', 14, 'Delimiter', ',');

tmp = file(13:end,8);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
temperaturePre = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)); % K. Remove missing elements
tmp = file(13:end,10);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
specHeatPre = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)).*molarMass./(mass./1e3)./(1e6); % Converting from muJ/K to J/mol-K
tmp = file(13:end,11);
mask = cellfun(@ismissing,tmp);
tmp(mask) = {[]};
specHeatErrPre = cell2mat(cellfun(@rmmissing, tmp, 'UniformOutput', false)).*molarMass./(mass./1e3)./(1e6); % Converting from muJ/K to J/mol-K

specHeatLBPre = specHeatPre*mass/(mass+massErr); % Lower envelope representing systematic uncertainty in mass.
specHeatLBErrPre = specHeatErrPre*mass/(mass+massErr);
specHeatUBPre = specHeatPre*mass/(mass-massErr); % Upper envelope representing systematic uncertainty in mass.
specHeatUBErrPre = specHeatErrPre*mass/(mass-massErr);

[specHeatCalc, specHeatMag, entropy, entropyErr, specHeat, specHeatErr, temperature] = anharmEntropy(specHeatPre, specHeatErrPre, temperaturePre, n, fitLimLow, fitLimHigh, debTemp0, anharm0, errPts, fact, offset);
%[specHeatLBCalc, specHeatLBMag, entropyLB, entropyLBErr, specHeatLB, specHeatLBErr, temperatureLB] = anharmEntropy(specHeatLBPre, specHeatLBErrPre, temperaturePre, n, fitLimLow, fitLimHigh, debTemp0, anharm0, errPts, fact, offset);
%[specHeatUBCalc, specHeatUBMag, entropyUB, entropyUBErr, specHeatUB, specHeatUBErr, temperatureUB] = anharmEntropy(specHeatUBPre, specHeatUBErrPre, temperaturePre, n, fitLimLow, fitLimHigh, debTemp0, anharm0, errPts, fact, offset);

figure(1)
nexttile(1)
hold on
%plot(temperatureLB, specHeatLB./temperatureLB, '--k') % Data.  Pretty much on top of the fit because the relative uncertainty in the mass is small.
%plot(temperatureUB, specHeatUB./temperatureUB, '--k')
xline(transition7K, '--k', 'LineWidth', 1.0); % Midpoint of TN2
xline(transition14K, '--k', 'LineWidth', 1.0); % Midpoint of TN1
xlim([0 270])
set(gca, 'ytick', [0, 2, 4, 6, 8, 10])
set(gca, 'xticklabel', [])
set(gca, 'XScale', 'log')
text(0.02, 0.9, '(a)', 'units', 'normalized', 'FontSize', 12)
hold off
nexttile(2)
hold on
%plot(temperatureLB, specHeatLBMag./temperatureLB, '--k') % Model. Pretty much on top of the fit because the relative uncertainty in the mass is small.
%plot(temperatureUB, specHeatUBMag./temperatureUB, '--k')
xline(transition7K, '--k', 'LineWidth', 1.0) % Midpoint of TN2
xline(transition14K, '--k', 'LineWidth', 1.0) % Midpoint of TN1
yline(0, '--k', 'LineWidth', 0.5) % x-axis to indicate which points are negative
ylim([-0.5 10])
xlim([0 270])
set(gca, 'ytick', [0, 2, 4, 6, 8])
set(gca, 'xticklabel', [])
set(gca, 'XScale', 'log')
text(0.02, 0.9, '(b)', 'units', 'normalized', 'FontSize', 12)
hold off
nexttile(3)
hold on
%plot(temperatureLB, entropyLB, '--k') % Model. Pretty much on top of the fit because the relative uncertainty in the mass is small.
%plot(temperatureUB, entropyUB, '--k')
ln = yline(5*N*kB*log(8),'--k', '\rm5\itR\rmln(8)', 'LabelHorizontalAlignment', 'center'); % 5 Eu atoms per formula unit, spin 7/2
ln.FontSize = 10;
xline(transition7K, '--k', 'LineWidth', 1.0) % Midpoint of TN2
xline(transition14K, '--k', 'LineWidth', 1.0) % Midpoint of TN1
xlim([0 270])
ylim([0 115])
set(gca, 'ytick', [0, 40, 80, 120])
set(gca, 'XScale', 'log')
text(0.02, 0.9, '(c)', 'units', 'normalized', 'FontSize', 12)
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\HC.eps';
%exportgraphics(gcf, fdir, 'ContentType', 'vector')

%save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\HC.mat') % Save the Workspace

% Estimate error from integral rather than cumtrapz. See 04/06/2023 526 notes.
entropyErrEst = sqrt(sum((temperature(2:end)-temperature(1:end-1)).^2./temperature(1:end-1).^2.*specHeatErr(1:end-1).^2));

%% Figure 3: Overplot
clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 10)

% Add the heat capacity plot
mass = 5.19; % Sample mass in mg
molarMass = 5*(151.96)+2*(114.82)+6*(121.76);

file = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Excel\Eu5In2Sb6\HC\HC12252020.xlsx');
temperature = rmmissing(file(1).data(:,8)); % K
specHeat = rmmissing(file(1).data(:,10))*molarMass/(mass/1e3)/(1e6); % Converting from muJ/K to J/mol-K
specHeatErr = rmmissing(file(1).data(:,11))*molarMass/(mass/1e3)/(1e6); % Converting from muJ/K to J/mol-K

% Average datapoints at essentially the same temperature
% nAvg = 3;
% temperature = arrayfun(@(i) mean(temperatureAll(i:i+nAvg-1)),1:nAvg:length(temperatureAll)-nAvg+1)';
% specHeat = arrayfun(@(i) mean(specHeatAll(i:i+nAvg-1)),1:nAvg:length(specHeatAll)-nAvg+1)';
% specHeatErr = arrayfun(@(i) sqrt(sum(specHeatErrAll(i:i+nAvg-1).^2))/length(specHeatErrAll(i:i+nAvg-1)),1:nAvg:length(specHeatErrAll)-nAvg+1)';

figure('Units', 'inches', 'Position', [0, 1.0, 3.375, 6.0])
t = tiledlayout(5, 1, 'TileSpacing', 'none', 'Padding', 'compact');
ha(1) = nexttile;
hold on
text(0.04, 0.9, '(a)', 'Units', 'normalized', 'fontsize', 11)
ylabel('\itC\rm (J/mol\cdotK)')
errorbar(temperature, specHeat, specHeatErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'XTickLabel', []);
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
box on
xlim([0 20])
ylim([0 140])
xline(7.2028, '--k', 'LineWidth', 1.0)
xline(14.240, '--k', 'LineWidth', 1.0)
xticks([0 5 10 15 20])
yticks([0 25 50 75 100 125])
hold off

% (1,0,0) OP
clear

% Fit for beta with chosen TMin from Eu5In2Sb6Beta.m result, TMin = 13.13 K
cutoff = 13.1; % Minimum temperature cutoffs to test
estTN = 14.23; % To determine cutoffs on data
fact = [0.5, 1, 0.01, 0];
offset = [0, 0, 0, 0.5];
x0 = [5, 0.5]; % Setting initial beta to mean-field so I don't bias the fitting too much.

% Import data as tables with variable names
filename1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\100_Tscan_854308.bt7';
filename2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\100_Tscan_854310.bt7';

data1 = readtable(filename1, 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 44);
fid1 = fopen(filename1); % Open file to take column headings as one line of text and edit
title1 = textscan(fid1, '%s', 145, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 43);
fclose(fid1);
title1 = title1{1};
data1.Properties.VariableNames = title1(2:end); % Remove "#Columns" from variable names

data2 = readtable(filename2, 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 44);
fid2 = fopen(filename2); % Open file to take column headings as one line of text and edit
title2 = textscan(fid2, '%s', 145, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 43);
fclose(fid2);
title2 = title2{1};
data2.Properties.VariableNames = title2(2:end); % Remove "#Columns" from variable names

data = [data1; data2]; % Concatenate tables

% Extract parameters
temp = data.Temp;
mon = data.Monitor;
time = data.Time;
det = data.Detector;
intPrelim = det./mon; % Normalize to monitor
intPrelimErr = intPrelim.*sqrt(1./det + 1./mon);
int = intPrelim.*mean(mon./time); % Convert to cts/sec given the typical mon/time ratio over the course of the scan
monAvgErr = sqrt(sum(mon./time))/length(mon); % The average of all of the monitor counts has some errorbar associated with it. Determined from formula for propagating uncorrelated errors.
intErr = int.*sqrt(intPrelimErr.^2./intPrelim.^2 + monAvgErr^2./mean(mon)^2); % Propagate errors from multiplying intPrelim with average of monitor counts then dividing by time. See OverplotErr in 080322 notes.

tempRed = @(TC, T) 1-T./TC; % Reduced temperature
model = @(x, T) x(1) + x(2).*tempRed(x(3), T).^(2.*x(4)).*(T<x(3)); % bg, I0, TC, beta. Flat background and critical scattering, see PAN manual from DAVE except for factor of 2.

tempFit = temp(temp>cutoff);
intFit = int(temp>cutoff);
intErrFit = intErr(temp>cutoff);

modelInput = @(x) model(x, tempFit);
x0 = [mean(intFit(tempFit>estTN)), x0(1), estTN, x0(2)];
errPts = 5e2;
[xFit, redChi2Fit, xErr, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(intFit, intErrFit, modelInput, x0, errPts, fact, offset);

tempCalc = linspace(tempFit(1), tempFit(end), 5e2); % Solid line over fitted region
tempCalcFull = linspace(7.3, tempFit(1), 5e2); % Dashed line beyond fitted region
intCalc = model(xFit(:), tempCalc);
intCalcFull = model(xFit(:), tempCalcFull);

disp(['Beta: ', num2str(xFit(4))])
disp(['Beta error: ', num2str(xErr(4))])

figure(1)
ha(2) = nexttile;
hold on
text(0.04, 0.1, '(b)', 'Units', 'normalized', 'fontsize', 11)
ylabel('\itI\rm (cts/sec.)')
errorbar(temp, int, intErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
p1 = plot(tempCalc, intCalc, 'LineWidth', 1, 'Color', 'r');
p2 = plot(tempCalcFull, intCalcFull, '--', 'LineWidth', 1, 'Color', 'r');
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'XTickLabel', []);
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
xlim([0 20])
ylim([0 (max(int)+max(intErr))*1.1])
xline(7.2028, '--k', 'LineWidth', 1.0)
xline(14.240, '--k', 'LineWidth', 1.0)
xticks([0 5 10 15 20])
yticks([0 5 10 15 20])
box on
text(0.95, 0.75, {'(100)', 'Peak', 'BT7'}, 'Units', 'normalized', 'fontsize', 9, 'HorizontalAlignment', 'right')
hold off

% (1,4,0) OP
clear

% Import data as tables with variable names
filename1 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\EUInS004.ng5';
filename2 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\EuInS047.ng5';
filename3 = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\EuInS048.ng5';

data1 = readtable(filename1, 'FileType', 'delimitedtext', 'NumHeaderLines', 11);
fid1 = fopen(filename1); % Open file to take column headings as one line of text and edit
data2 = readtable(filename2, 'FileType', 'delimitedtext', 'NumHeaderLines', 11);
fid2 = fopen(filename2); % Open file to take column headings as one line of text and edit
data3 = readtable(filename3, 'FileType', 'delimitedtext', 'NumHeaderLines', 11);
fid3 = fopen(filename3); % Open file to take column headings as one line of text and edit

data = [data1; data2; data3]; % Concatenate tables

% Extract parameters
temp = data.T_act;
mon = 4.67e6; % Monitor counts for each datapoint. Can confirm in DAVE the mon has to be hit by prf.
time = data.min*60;
det = data.Counts;
intPrelim = det./mon; % Normalize to monitor
intPrelimErr = intPrelim.*sqrt(1./det + 1./mon);
int = intPrelim.*mean(mon./time); % Convert to cts/sec given the typical mon/time ratio over the course of the scan
monAvgErr = sqrt(sum(mon./time))/length(mon); % The average of all of the monitor counts has some errorbar associated with it. Determined from formula for propagating uncorrelated errors.
intErr = int.*sqrt(intPrelimErr.^2./intPrelim.^2 + monAvgErr^2./mean(mon)^2); % Propagate errors from multiplying intPrelim with average of monitor counts then dividing by time. See OverplotErr in 080322 notes.

ha(3) = nexttile;
hold on
text(0.04, 0.1, '(c)', 'Units', 'normalized', 'fontsize', 11)
ylabel('\itI\rm (cts/sec.)')
errorbar(temp, int, intErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
set(gca, 'TickLength', [0.02, 0.01])
xlim([0 20])
ylim([4 8])
xline(7.2028, '--k', 'LineWidth', 1.0)
xline(14.240, '--k', 'LineWidth', 1.0)
set(gca,'XTickLabel',[]);
set(gca,'PlotBoxAspectRatio',[2 1 1]);
xticks([0 5 10 15 20])
yticks([4 5 6 7])
box on
text(0.95, 0.75, {'(140)', 'Peak', 'SPINS'}, 'Units', 'normalized', 'fontsize', 9, 'HorizontalAlignment', 'right')
hold off

% (0,1,1) OP
clear

filename = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\011_Tscan_854230.bt7';

data = readtable(filename, 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 44);
fid = fopen(filename); % Open file to take column headings as one line of text and edit
heading = textscan(fid, '%s', 145, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 43);
fclose(fid);
heading = heading{1};
data.Properties.VariableNames = heading(2:end); % Remove "#Columns" from variable names
data(37:end,:) = []; % Delete rows after neutron source failed

% Extract parameters
temp = data.Temp;
mon = data.Monitor;
time = data.Time;
det = data.Detector;
intPrelim = det./mon; % Normalize to monitor
intPrelimErr = intPrelim.*sqrt(1./det + 1./mon);
int = intPrelim.*mean(mon./time); % Convert to cts/sec given the typical mon/time ratio over the course of the scan
monAvgErr = sqrt(sum(mon./time))/length(mon); % The average of all of the monitor counts has some errorbar associated with it. Determined from formula for propagating uncorrelated errors.
intErr = int.*sqrt(intPrelimErr.^2./intPrelim.^2 + monAvgErr^2./mean(mon)^2); % Propagate errors from multiplying intPrelim with average of monitor counts then dividing by time. See OverplotErr in 080322 notes.

ha(4) = nexttile;
hold on
text(0.04, 0.1, '(d)', 'Units', 'normalized', 'fontsize', 11)
ylabel('\itI\rm (cts/sec.)')
errorbar(temp, int, intErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
set(gca, 'TickLength', [0.02, 0.01])
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
xlim([0 20])
ylim([0 (max(int)+max(intErr))*1.1])
xline(7.2028, '--k', 'LineWidth', 1.0)
xline(14.240, '--k', 'LineWidth', 1.0)
set(gca,'XTickLabel',[]);
xticks([0 5 10 15 20])
yticks([0 1 2])
box on
text(0.95, 0.75, {'(011)', 'Peak', 'BT7'}, 'Units', 'normalized', 'fontsize', 9, 'HorizontalAlignment', 'right')
hold off

% (0,0,1/2) OP
clear

errPtsSig = 1e2; % For a3/a4 sigma, center fits
factSig = [0, 0.3, 0.3, 0]; % bg, area, sigma, center
offsetSig = [0.01, 0.04, 1e-3, 0.1];
errPtsPrelim = 1e3; % For a3/a4 area fits
factPrelim = [1, 1]; % bg, area
offsetPrelim = [0, 1];

% Import all a3/a4 scans as a function of temperature
filenums = 854258:854303;
filenames = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\fpx'+string(filenums)+'.bt7';
for i = 1:length(filenums)
    data{i} = readtable(filenames(i), 'FileType', 'delimitedtext', 'TreatAsMissing', 'N/A', 'NumHeaderLines', 44);
    fid = fopen(filenames(i)); % Open file to take column headings as one line of text and edit
    heading = textscan(fid, '%s', 145, 'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'headerlines', 43);
    fclose(fid);
    heading = heading{1};
    data{i}.Properties.VariableNames = heading(2:end); % Remove "#Columns" from variable names
    data{i}(end,:) = [];

    temp(i) = mean(data{i}.Temp);
    a4(:,i) = data{i}.A4;
    mon(:,i) = data{i}.Monitor;
    time(:,i) = data{i}.Time;
    det(:,i) = data{i}.Detector;
    intPrelim(:,i) = det(:,i)./mon(:,i); % Normalize to monitor
    intPrelimErr(:,i) = intPrelim(:,i).*sqrt(1./det(:,i) + 1./mon(:,i));
    int(:,i) = intPrelim(:,i).*mean(mon(:,i)./time(:,i)); % Convert to cts/sec given the typical mon/time ratio over the course of the scan
    monAvgErr(:,i) = sqrt(sum(mon(:,i)./time(:,i)))/length(mon(:,i)); % The average of all of the monitor counts has some errorbar associated with it. Determined from formula for propagating uncorrelated errors.
    intErr(:,i) = int(:,i).*sqrt(intPrelimErr(:,i).^2./intPrelim(:,i).^2 + monAvgErr(:,i)^2./mean(mon(:,i))^2); % Propagate errors from multiplying intPrelim with average of monitor counts then dividing by time. See OverplotErr in 080322 notes.
end

%First fix sigma based on average of first 5 peaks to better-control the
%high-T fits. Gaussian and flat background.
modelSig = @(x, a4) x(1) + x(2)./sqrt(2.*pi)./x(3).*exp(-((a4-x(4))./x(3)).^2./2); % bg, area, sigma, center
for i = 1:5
    modelInputSig = @(x) modelSig(x, a4(:,i));
    tmp = sort(int(:,i)); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(a4(:,i).*(int(:,i)-bg0))./sum(int(:,i)-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    sig0 = 0.2;
    area0 = (max(int(:,i))-bg0)*sqrt(2*pi)*sig0; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0Sig = [bg0, area0, sig0, cent0];
    [xFitSig(:,i), redChi2FitSig(:,i), xErrSig(:,i), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(int(:,i), intErr(:,i), modelInputSig, x0Sig, errPtsSig, factSig, offsetSig);

    % Plot a4 fit
    a4Cal = linspace(min(a4(:,i)), max(a4(:,i)), 5e2)';
    intCal = modelSig(xFitSig(:,i), a4Cal);
    figure(2)
    close
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title('(0 0 1/2) A_4 Scan: \sigma')
    xlabel('A_4 (deg.)')
    ylabel('Intensity (det. deg./mon.)')
    errorbar(a4(:,i), int(:,i), intErr(:,i), 'o', 'MarkerFaceColor', 'w')
    plot(a4Cal, intCal)
    box on
    axis square
    hold off
    pause(0.1)
end
sig = mean(xFitSig(3,:));
center = mean(xFitSig(4,:));

% Fix the width and fit for integrated intensity. Gaussian and flat
% background
modelPrelim = @(x, a4) x(1) + x(2)./sqrt(2.*pi)./sig.*exp(-((a4-center)./sig).^2./2); % bg, area
for i = 1:length(filenums)
    modelInputPrelim = @(x) modelPrelim(x, a4(:,i));
    tmp = sort(int(:,i)); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    area0 = (max(int(:,i))-bg0)*sqrt(2*pi)*sig; % Guess that the area is given by its relation in terms of the assumed peak height and assumed standard deviation where the peak height is estimated as the largest observed intensity in the scan minus the guessed background.
    x0Prelim = [bg0, area0];
    [xFitPrelim(:,i), redChi2FitPrelim(:,i), xErrPrelim(:,i), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(int(:,i), intErr(:,i), modelInputPrelim, x0Prelim, errPtsPrelim, factPrelim, offsetPrelim);
    area(:,i) = xFitPrelim(2,i); % Integrated intensity
    areaErr(:,i) = xErrPrelim(2,i);

    % Plot a4 fit
    a4Cal = linspace(min(a4(:,i)), max(a4(:,i)), 5e2)';
    intCal = modelPrelim(xFitPrelim(:,i), a4Cal);
    figure(2)
    close
    figure('Units', 'normalized', 'Position', [0, 0.3, 0.5, 0.6])
    clf
    hold on
    title('(0 0 1/2) A_4 Scan')
    xlabel('A_4 (deg.)')
    ylabel('Intensity')
    errorbar(a4(:,i), int(:,i), intErr(:,i), 'o', 'MarkerFaceColor', 'w')
    plot(a4Cal, intCal)
    box on
    axis square
    hold off
    pause(0.1)
end
figure(2)
close

% Calculate expected critical behavior for 3D Heisenbergy beta=0.366, Tc=7.2028
cutoff = 5.8;
tempFit = temp(temp>cutoff);
intFit = area(temp>cutoff);
intErrFit = areaErr(temp>cutoff);

tempRed = @(TC, T) 1-T./TC; % Reduced temperature
model = @(x, Tc, beta, T) x(1) + x(2).*tempRed(Tc, T).^(2.*beta).*(T<Tc); % bg, I0. Flat background and critical scattering, see PAN manual from DAVE except for factor of 2.

fact = [1, 1];
offset = [0, 1];
modelInput = @(x) model(x, 7.2028, 0.366, tempFit);
x0 = [mean(intFit(tempFit>7.2028)), 1];
errPts = 5e2;
[xFit, redChi2Fit, xErr, ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(intFit, intErrFit, modelInput, x0, errPts, fact, offset);

tempCalc = linspace(tempFit(1), tempFit(end), 5e2); % Solid line over fitted region
tempCalcFull = linspace(4.0, tempFit(1), 5e2); % Dashed line beyond fitted region
intCalc = model(xFit(:), 7.2028, 0.366, tempCalc);
intCalcFull = model(xFit(:), 7.2028, 0.366, tempCalcFull);

ha(5) = nexttile;
hold on
text(0.04, 0.1, '(e)', 'Units', 'normalized', 'fontsize', 11)
ylabel('\itI\rm (cts\cdotdeg.(2\Theta)/sec.)')
xlabel('\itT\rm (K)')
errorbar(temp, area, areaErr, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
p1 = plot(tempCalc, intCalc, 'LineWidth', 1, 'Color', 'r');
p2 = plot(tempCalcFull, intCalcFull, '--', 'LineWidth', 1, 'Color', 'r');
set(gca, 'TickLength', [0.02, 0.01])
xlim([0 20])
ylim([0 (max(area)+max(areaErr))*1.1])
xline(7.2028, '--k', 'LineWidth', 1.0)
xline(14.240, '--k', 'LineWidth', 1.0)
set(gca,'PlotBoxAspectRatio',[2 1 1]);
xticks([0 5 10 15 20])
yticks([0 5 10 15 20])
box on
text(0.95, 0.75, {'(0,0,1/2)', '\Theta-2\Theta', 'BT7'}, 'Units', 'normalized', 'fontsize', 9, 'HorizontalAlignment', 'right')
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\overplot.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%% Figure 4: Peaks
clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 10)

% Fit to the (100) reflection from SPINS experiment with Voigt and pick out
% Lorentzian hwhm for correlation length

% Variables for correlation length fit
sigSPINS = 0.48742/(2*sqrt(2*log(2))); % See Neutron Analysis script
a = 12.51000; %A: angstroms ICSD_CollCode281238.cif
b = 14.58400; %B: angstroms
c = 4.62430; %C: angstroms
HKL = [1, 0, 0];
q = 2*pi*sqrt(HKL(1)^2/a^2 + HKL(2)^2/b^2 + HKL(3)^2/c^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
errPts = 0; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.3, 0.005, 0.0]; % Determines range of iterated values for errorbar. [bg, peak, center, gamma]
offset = [0.0, 0.2, 0.005, 0.1]; % Determines range of iterated values for errorbar
model = @(x, a3, sig) x(1) + x(2).*voigt(a3, x(3), sig, x(4)); % Voigt function with flat background.

% To convert to reciprocal lattice units later
aStar = 2.*pi./a;
bStar = 2.*pi./b;
cStar = 2.*pi./c;

% Import data
file1 = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\fpx03063.ng5');
file3 = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\fpx03064.ng5');
file2 = importdata('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\092920Data\fpx03065.ng5');

motor = [file1.data(:,1), file2.data(:,1), file3.data(:,1)];
int = [file1.data(:,2), file2.data(:,2), file3.data(:,2)];
intErr = sqrt(int);

% Fit for correlation length
for i = 1:length(motor(1,:))
    tmp = sort(int(:,i)); % Sorted intensities
    bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(motor(:,i).*(int(:,i)-bg0))./sum(int(:,i)-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    gam0 = 0.01;
    peak0 = (max(int(:,i))-bg0)/voigt(cent0, cent0, sigSPINS, gam0); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
    x0 = [bg0, peak0, cent0, gam0];
    modelInput = @(x) model(x, motor(:,i), sigSPINS);
    [xFit(:,i), redChi2Fit(:,i), xErr(:,i), chiUpper(:,i), chiTrial(:,:,i), paramTrial(:,:,i), interpPts(:,:,i), ~, ~, paramLower(:,i), paramUpper(:,i)] = fitRedChi2Err(int(:,i), intErr(:,i), modelInput, x0, errPts, fact, offset);
    fwhm(i) = 2*q*xFit(4,i)*pi/180; % Fwhm from Lorentzian component of Voigt function
    fwhmErr(i) = 2*q*xErr(4,i)*pi/180;
    corr(i) = 2/fwhm(i); % Factor of 2 is important to remember
    corrErr(i) = 2*fwhmErr(i)/fwhm(i)^2; % Square is from the quadrature sum partial derivatives times errors
    fwhmMin(i) = 2*q*paramUpper(4, i)*pi/180; % Using fwhm max b/c refers to corr min
    corrMin(i) = 2/fwhmMin(i); % Factor of 2 is important to remember
    fwhmMax(i) = 2*q*paramLower(4, i)*pi/180; % Fwhm from Lorentzian component of Voigt function
    corrMax(i) = 2/fwhmMax(i); % Factor of 2 is important to remember

    % Values for plotting the fit
    motorCal(:,i) = linspace(min(motor(:,i)), max(motor(:,i)), 500)';
    intCal(:,i) = model(xFit(:,i), motorCal(:,i), sigSPINS);

    % Convert to rlu, in HK0 looking at 100 with a3 scan so taking 0K0 projection
    qK(:,i) = q.*sind(mean(motor(:,i), 1).*ones(size(motor(:,i)))-motor(:,i)); % Here the a3 angle in real space doesn't correspond to the rotated angle in reciprocal space because the lattice is orthorhombic.
    K(:,i) = qK(:,i)./bStar;
    qKCal(:,i) = q.*sind(mean(motorCal(:,i), 1).*ones(size(motorCal(:,i)))-motorCal(:,i));
    KCal(:,i) = qKCal(:,i)./bStar;
end

% Plot scans and fit, display correlation length, display errorbars in a
% separate figure for reference
figure('Units', 'inches', 'Position', [0, 1.0, 3.375, 3.78])
tiledlayout(2, 1, 'TileSpacing', 'compact')
nexttile
hold on
e1=errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
e2=errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
e3=errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%p1=plot(KCal(:,1), intCal(:,1), 'Color', 'r');
p2=plot(KCal(:,2), intCal(:,2), 'Color', 'b', 'LineWidth', 1);
p3=plot(KCal(:,3), intCal(:,3), 'Color', 'k', 'LineWidth', 1);
legend([e1, e2, e3],{'\itT\rm = 20 K', '\itT\rm = 10 K', '\itT\rm = 1.5 K'})
legend('boxoff')
text(0.15, 0.9, '(100)', 'Units', 'normalized', 'fontsize', 11)
text(0.15, 0.75, 'SPINS', 'Units', 'normalized', 'fontsize', 11)
%xlabel('(0K0) (r.l.u.)')
ylabel('\itI\rm (det. cts/7e5 mon. cts)')
xlim([min(K, [], 'all'), max(K, [], 'all')])
ylim([0 2200])
yticks(linspace(0, 1700, 5))
set(gca, 'TickLength',[0.02, 0.025]) % [2Dlength 3Dlength]
pbaspect([16 9 1])
box on
text(0.02, 0.9, '(a)', 'units', 'normalized', 'fontsize', 12)
hold off

%disp(['20 K (1,0,0) Correlation Length: ', num2str(corr(1)), '+-', num2str(corrErr(1)), ' Angstroms'])
disp(['10 K (1,0,0) Correlation Length: ', num2str(corr(2)), '+-', num2str(corrErr(2)), ' Angstroms'])
disp(['10 K (1,0,0) Correlation Length Min: ', num2str(corrMin(2)), ' Angstroms'])
disp(['10 K (1,0,0) Correlation Length Max: ', num2str(corrMax(2)), ' Angstroms'])
disp(['1.5 K (1,0,0) Correlation Length: ', num2str(corr(3)), '+-', num2str(corrErr(3)), ' Angstroms'])
disp(['1.5 K (1,0,0) Correlation Length Min: ', num2str(corrMin(3)), ' Angstroms'])
disp(['1.5 K (1,0,0) Correlation Length Max: ', num2str(corrMax(3)), ' Angstroms'])

tempText = {'20 K', '10 K', '1.5 K'};
for i = 2:length(motor(1,:))
    figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
    for j = 1:length(x0)
        subplot(2, 2, j)
        hold on
        ylabel('\chi^2_{r}')
        if ~any(isnan(interpPts(:, j, i))) % Check to see if an errorbar was determined before plotting
            scatter(paramTrial(:, j, i), chiTrial(:, j, i))
            plot([paramTrial(interpPts(1, j, i), j, i), paramTrial(interpPts(2, j, i), j, i)], [chiTrial(interpPts(1, j, i), j, i), chiTrial(interpPts(2, j, i), j, i)], 'b')
            plot([paramTrial(interpPts(3, j, i), j, i), paramTrial(interpPts(4, j, i), j, i)], [chiTrial(interpPts(3, j, i), j, i), chiTrial(interpPts(4, j, i), j, i)], 'b')
            yline(chiUpper, 'Color', 'r', 'LineWidth', 3.0)
            plot([paramLower(j, i), paramUpper(j, i)], [chiUpper(i), chiUpper(i)], 'k-.o', 'LineWidth', 2.0)
            xlim([-inf inf])
            ylim([redChi2Fit(i)-0.1 chiUpper(i)+(chiUpper(i)-redChi2Fit(i))])
        end
        set(gca,'FontSize',12)
        box on
        hold off
    end
    subplot(2, 2, 1)
    hold on
    title([tempText{i}, ' (100) Background: ', num2str(round(xFit(1, i), 4)), '\pm', num2str(round(xErr(1, i), 4))])
    xlabel('Background (det. cts/7e5 mon. cts)')
    hold off
    subplot(2, 2, 2)
    hold on
    title([tempText{i}, ' (100) Voigt Peak: ', num2str(round(xFit(2, i), 4)), '\pm', num2str(round(xErr(2, i), 4))])
    xlabel('Voigt Peak (det. cts/7e5 mon. cts)')
    hold off
    subplot(2, 2, 3)
    hold on
    title([tempText{i}, ' (100) Voigt Center: ', num2str(round(xFit(3, i), 4)), '\pm', num2str(round(xErr(3, i), 4))])
    xlabel('Voigt Center (deg.)')
    hold off
    subplot(2, 2, 4)
    hold on
    title([tempText{i}, ' (100) Lorentzian FWHM: ', num2str(round(xFit(4, i), 5)), '\pm', num2str(round(xErr(4, i), 5))])
    xlabel('Lorentzian FWHM (deg.)')
    hold off
end

%save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\CorrBestSPINS.mat') % Save the Workspace

% Add the (0,0,3/2) reflection from BT7 experiment
% Fit to the (0,0,3/2) reflection from SPINS experiment with Voigt and pick out
% Lorentzian hwhm for correlation length
clear

% Variables for correlation length fit
sigBT7 = 0.38748/(2*sqrt(2*log(2)));
a = 12.51000; %A: angstroms ICSD_CollCode281238.cif
b = 14.58400; %B: angstroms
c = 4.62430; %C: angstroms
HKL = [0, 0, 3/2];
q = 2*pi*sqrt(HKL(1)^2/a^2 + HKL(2)^2/b^2 + HKL(3)^2/c^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
errPts = 0; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.3, 0.005, 0.0]; % Determines range of iterated values for errorbar. [bg, peak, center, gamma]
offset = [0.0, 0.2, 0.005, 1e-8]; % Determines range of iterated values for errorbar
model = @(x, a3, sig) x(1) + x(2).*voigt(a3, x(3), sig, x(4)); % Voigt function with flat background.

% To convert to reciprocal lattice units later
aStar = 2.*pi./a;
bStar = 2.*pi./b;
cStar = 2.*pi./c;

file1 = readcell(strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\fpx854193.bt7'), 'FileType', 'text', 'NumHeaderLines', 44, 'CommentStyle', '#', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
file2 = readcell(strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\fpx854147.bt7'), 'FileType', 'text', 'NumHeaderLines', 44, 'CommentStyle', '#', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');

motor = [cell2mat(file1(:,13)), cell2mat(file2(:,13))];
int = [cell2mat(file1(:,10)), cell2mat(file2(:,10))];
intErr = sqrt(int);

% Fit for correlation length
for i = 1:length(motor(1,:))
    tmp = sort(int(:,i)); % Sorted intensities
    bg0 = mean(tmp(1:6)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
    cent0 = sum(motor(:,i).*(int(:,i)-bg0))./sum(int(:,i)-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
    gam0 = 0.01;
    peak0 = (max(int(:,i))-bg0)/voigt(cent0, cent0, sigBT7, gam0); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
    x0 = [bg0, peak0, cent0, gam0];
    modelInput = @(x) model(x, motor(:,i), sigBT7);
    [xFit(:,i), redChi2Fit(:,i), xErr(:,i), chiUpper(:,i), chiTrial(:,:,i), paramTrial(:,:,i), interpPts(:,:,i), ~, ~, paramLower(:,i), paramUpper(:,i)] = fitRedChi2Err(int(:,i), intErr(:,i), modelInput, x0, errPts, fact, offset);
    fwhm(i) = 2*q*xFit(4,i)*pi/180; % Fwhm from Lorentzian component of Voigt function
    fwhmErr(i) = 2*q*xErr(4,i)*pi/180;
    corr(i) = 2/fwhm(i); % Factor of 2 is important to remember
    corrErr(i) = 2*fwhmErr(i)/fwhm(i)^2; % Square is from the quadrature sum partial derivatives times errors
    fwhmMin(i) = 2*q*paramUpper(4, i)*pi/180; % Using fwhm max b/c refers to corr min
    corrMin(i) = 2/fwhmMin(i); % Factor of 2 is important to remember
    fwhmMax(i) = 2*q*paramLower(4, i)*pi/180; % Fwhm from Lorentzian component of Voigt function
    corrMax(i) = 2/fwhmMax(i); % Factor of 2 is important to remember

    motorCal(:,i) = linspace(min(motor(:,i)), max(motor(:,i)), 500)';
    intCal(:,i) = model(xFit(:,i), motorCal(:,i), sigBT7);

    % Convert to rlu, in 0KL looking at 001p5 with a3 scan so taking 0K0 projection
    qK(:,i) = q.*sind(mean(motor(:,i), 1).*ones(size(motor(:,i)))-motor(:,i)); % Here the a3 angle in real space doesn't correspond to the rotated angle in reciprocal space because the lattice is orthorhombic.
    K(:,i) = qK(:,i)./bStar;
    qKCal(:,i) = q.*sind(mean(motorCal(:,i), 1).*ones(size(motorCal(:,i)))-motorCal(:,i));
    KCal(:,i) = qKCal(:,i)./bStar;
end

% Plot scans and fit, display correlation length, display errorbars in a
% separate figure for reference
figure(1)
nexttile
hold on
e1 = errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
e2 = errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%p1 = plot(KCal(:,1), intCal(:,1), 'Color', 'r');
p2 = plot(KCal(:,2), intCal(:,2), 'Color', 'k', 'LineWidth', 1);
legend([e1, e2],{'\itT\rm = 18 K', '\itT\rm = 1.6 K'})
legend('boxoff')
text(0.15, 0.9, '(0,0,3/2)', 'Units', 'normalized', 'fontsize', 11)
text(0.15, 0.75, 'BT7', 'Units', 'normalized', 'fontsize', 11)
xlabel('(0\itK\rm0) (rlu)')
ylabel('\itI\rm (det. cts/15 sec.)')
xlim([min(K, [], 'all'), max(K, [], 'all')])
ylim([0 1600])
yticks(linspace(0, 1400, 5))
set(gca, 'TickLength',[0.02, 0.025]) % [2Dlength 3Dlength]
pbaspect([16 9 1])
box on
text(0.02, 0.9, '(b)', 'units', 'normalized', 'fontsize', 12)
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\peaks.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%disp(['18 K (0,0,3/2) Correlation Length: ', num2str(corr(1)), '+-', num2str(corrErr(1)), ' Angstroms'])
disp(['1.6 K (0,0,3/2) Correlation Length: ', num2str(corr(2)), '+-', num2str(corrErr(2)), ' Angstroms'])
disp(['1.6 K (0,0,3/2) Correlation Length Min: ', num2str(corrMin(2)), ' Angstroms'])
disp(['1.6 K (0,0,3/2) Correlation Length Max: ', num2str(corrMax(2)), ' Angstroms'])

tempText = {'18 K', '1.6 K'};
i = 2;
figure('Units', 'normalized', 'Position', [0.5, 0.3, 0.5, 0.6])
for j = 1:length(x0)
    subplot(2, 2, j)
    hold on
    ylabel('\chi^2_{r}')
    if ~any(isnan(interpPts(:, j, i))) % Check to see if an errorbar was determined before plotting
        scatter(paramTrial(:, j, i), chiTrial(:, j, i))
        plot([paramTrial(interpPts(1, j, i), j, i), paramTrial(interpPts(2, j, i), j, i)], [chiTrial(interpPts(1, j, i), j, i), chiTrial(interpPts(2, j, i), j, i)], 'b')
        plot([paramTrial(interpPts(3, j, i), j, i), paramTrial(interpPts(4, j, i), j, i)], [chiTrial(interpPts(3, j, i), j, i), chiTrial(interpPts(4, j, i), j, i)], 'b')
        yline(chiUpper, 'Color', 'r', 'LineWidth', 3.0)
        plot([paramLower(j, i), paramUpper(j, i)], [chiUpper(i), chiUpper(i)], 'k-.o', 'LineWidth', 2.0)
        xlim([-inf inf])
        ylim([redChi2Fit(i)-0.1 chiUpper(i)+(chiUpper(i)-redChi2Fit(i))])
    end
    set(gca,'FontSize',12)
    box on
    hold off
end
subplot(2, 2, 1)
hold on
title([tempText{i}, ' (0,0,3/2) Background: ', num2str(round(xFit(1, i), 4)), '\pm', num2str(round(xErr(1, i), 4))])
xlabel('Background (det. cts/7e5 mon. cts)')
hold off
subplot(2, 2, 2)
hold on
title([tempText{i}, ' (0,0,3/2) Voigt Peak: ', num2str(round(xFit(2, i), 4)), '\pm', num2str(round(xErr(2, i), 4))])
xlabel('Voigt Peak (det. cts/7e5 mon. cts)')
hold off
subplot(2, 2, 3)
hold on
title([tempText{i}, ' (0,0,3/2) Voigt Center: ', num2str(round(xFit(3, i), 4)), '\pm', num2str(round(xErr(3, i), 4))])
xlabel('Voigt Center (deg.)')
hold off
subplot(2, 2, 4)
hold on
title([tempText{i}, ' (0,0,3/2) Lorentzian FWHM: ', num2str(round(xFit(4, i), 5)), '\pm', num2str(round(xErr(4, i), 5))])
xlabel('Lorentzian FWHM (deg.)')
hold off

%save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\CorrBestBT7.mat') % Save the Workspace

%% Figure 5: Refinements
clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')
set(0, 'defaultaxesfontsize', 14)

% Scale factors from paramagnetic refinements
scale20KHK0 = 2.210;
scale18K0KL = 0.6878;
scale10KHK0 = 2.210;
scale1p6K0KL = 0.6878;
scale1p5KHK0 = 2.210;
scale1p6K0KLHet = 0.6878;
scale1p5KHK0Het = 2.210;

% Import data
file20KHK0 = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\20KHK0\Eu5In2Sb620KHK0.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);
file18K0KL = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\18K0KL\Eu5In2Sb618K0KL.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);
file10KHK0 = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\10KHK0000Irrep3NoDomSubtDataB\Eu5In2Sb610KHK0Subt.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);
file1p6K0KL = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\Homogeneous\1p6k0KL000p5Irrep3NoDomAWin\PPP\Eu5In2Sb61p6K0KL000p5.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);
file1p5KHK0 = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\Homogeneous\1p5kHK0000Irrep3NoDomSubtDataB\Eu5In2Sb61p5KHK0Subt.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);
file1p6K0KLHet = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\Heterogeneous\1p6k0KL000p5Irrep3NoDomA\PPP\Eu5In2Sb61p6K0KL000p5.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);
file1p5KHK0Het = readcell('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\FullProf\Eu5In2Sb6\Finished\SphericalCorrection\Heterogeneous\1p5kHK0000Irrep3NoDomSubtDataB\PNN\Eu5In2Sb61p5KHK0Subt.prf', 'FileType', 'text', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join', 'NumHeaderLines', 2);

%Convert structure factors to barns
F2obs20KHK0 = cell2mat(file20KHK0(:,3))/scale20KHK0;
SigmaObs20KHK0 = cell2mat(file20KHK0(:,5))/scale20KHK0;
F2cal20KHK0 = cell2mat(file20KHK0(:,4))/scale20KHK0;
F2obs18K0KL = cell2mat(file18K0KL(:,3))/scale18K0KL;
SigmaObs18K0KL = cell2mat(file18K0KL(:,5))/scale18K0KL;
F2cal18K0KL = cell2mat(file18K0KL(:,4))/scale18K0KL;
F2obs10KHK0 = cell2mat(file10KHK0(:,3))/scale10KHK0;
SigmaObs10KHK0 = cell2mat(file10KHK0(:,5))/scale10KHK0;
F2cal10KHK0 = cell2mat(file10KHK0(:,4))/scale10KHK0;
F2obs1p6K0KL = cell2mat(file1p6K0KL(:,3))/scale1p6K0KL;
SigmaObs1p6K0KL = cell2mat(file1p6K0KL(:,5))/scale1p6K0KL;
F2cal1p6K0KL = cell2mat(file1p6K0KL(:,4))/scale1p6K0KL;
F2obs1p5KHK0 = cell2mat(file1p5KHK0(:,3))/scale1p5KHK0;
SigmaObs1p5KHK0 = cell2mat(file1p5KHK0(:,5))/scale1p5KHK0;
F2cal1p5KHK0 = cell2mat(file1p5KHK0(:,4))/scale1p5KHK0;
F2obs1p6K0KLHet = cell2mat(file1p6K0KLHet(:,3))/scale1p6K0KLHet;
SigmaObs1p6K0KLHet = cell2mat(file1p6K0KLHet(:,5))/scale1p6K0KLHet;
F2cal1p6K0KLHet = cell2mat(file1p6K0KLHet(:,4))/scale1p6K0KLHet;
F2obs1p5KHK0Het = cell2mat(file1p5KHK0Het(:,3))/scale1p5KHK0Het;
SigmaObs1p5KHK0Het = cell2mat(file1p5KHK0Het(:,5))/scale1p5KHK0Het;
F2cal1p5KHK0Het = cell2mat(file1p5KHK0Het(:,4))/scale1p5KHK0Het;

% Plot data
% Determine the lines with slope 1 to overplot. Make sure non-negative for
% log-log plots. Also make sure includes errorbars.
LinPM = linspace(min(nonzeros(abs([[F2cal20KHK0; F2cal18K0KL]; [F2obs20KHK0; F2obs18K0KL]-[SigmaObs20KHK0; SigmaObs18K0KL]])))*0.95, max([F2obs20KHK0; F2obs18K0KL]+[SigmaObs20KHK0; SigmaObs18K0KL])*1.05*100);
LinPh1 = linspace(min(nonzeros(abs([F2cal10KHK0; F2obs10KHK0-SigmaObs10KHK0])))*0.95, max([F2cal10KHK0; F2obs10KHK0+SigmaObs10KHK0])*1.05*100);
LinPh2 = linspace(min(nonzeros(abs([[F2cal1p6K0KL; F2cal1p5KHK0]; [F2obs1p6K0KL; F2obs1p5KHK0]-[SigmaObs1p6K0KL; SigmaObs1p5KHK0]])))*0.95, max([F2obs1p6K0KL; F2obs1p5KHK0]+[SigmaObs1p6K0KL; SigmaObs1p5KHK0])*1.05*100);
LinPhHet = linspace(min(nonzeros(abs([[F2cal1p6K0KLHet; F2cal1p5KHK0Het]; [F2obs1p6K0KLHet; F2obs1p5KHK0Het]-[SigmaObs1p6K0KLHet; SigmaObs1p5KHK0Het]])))*0.95, max([F2obs1p6K0KLHet; F2obs1p5KHK0Het]+[SigmaObs1p6K0KLHet; SigmaObs1p5KHK0Het])*1.05*100);

figure('Units', 'inches', 'Position', [0, 1.0, 7.2, 6.9])
t = tiledlayout(2, 2, 'TileSpacing', 'Compact');
nexttile
hold on
xlabel('\it{F}\rm{^2_{cal}(b)}')
ylabel('\it{F}\rm{^2_{obs}(b)}')
e20KHK0 = errorbar(F2cal20KHK0, F2obs20KHK0, SigmaObs20KHK0, 'o', 'MarkerFaceColor', 'w');
e18K0KL = errorbar(F2cal18K0KL, F2obs18K0KL, SigmaObs18K0KL, 'o', 'MarkerFaceColor', 'w');
plot(LinPM, LinPM, 'k', 'LineWidth', 1)
legend([e20KHK0, e18K0KL], {'20 K, (HK0)', '18 K, (0KL)'});
axis([min(LinPM), max(LinPM), min(LinPM), max(LinPM)])
axis square
set(gca,'XScale','log','YScale','log')
set(gca, 'XTick', [10^-2, 10^0, 10^2, 10^4])
set(gca, 'YTick', [10^-2, 10^0, 10^2, 10^4])
set(gca, 'TickLength', [0.04, 0.01])
box on
hold off

nexttile
hold on
xlabel('\it{F}\rm{^2_{cal}(b)}')
ylabel('\it{F}\rm{^2_{obs}(b)}')
e10KHK0 = errorbar(F2cal10KHK0, F2obs10KHK0, SigmaObs10KHK0, 'o', 'MarkerFaceColor', 'w');
plot(LinPh1, LinPh1, 'k', 'LineWidth', 1)
legend(e10KHK0, {'10 K, (HK0)'});
axis([min(LinPh1), max(LinPh1), min(LinPh1), max(LinPh1)])
axis square
set(gca,'XScale','log','YScale','log')
set(gca, 'XTick', [10^0, 10^2, 10^4])
set(gca, 'YTick', [10^0, 10^2, 10^4])
set(gca, 'TickLength', [0.04, 0.01])
box on
hold off

nexttile
hold on
xlabel('\it{F}\rm{^2_{cal}(b)}')
ylabel('\it{F}\rm{^2_{obs}(b)}')
e1p5KHK0 = errorbar(F2cal1p5KHK0, F2obs1p5KHK0, SigmaObs1p5KHK0, 'o', 'MarkerFaceColor', 'w');
e1p6K0KL = errorbar(F2cal1p6K0KL, F2obs1p6K0KL, SigmaObs1p6K0KL, 'o', 'MarkerFaceColor', 'w');
plot(LinPh2, LinPh2, 'k', 'LineWidth', 1)
legend([e1p5KHK0, e1p6K0KL], {'1.5 K, (HK0) Hom.', '1.6 K, (0KL) Hom.'});
axis([min(LinPh2), max(LinPh2), min(LinPh2), max(LinPh2)])
axis square
set(gca,'XScale','log','YScale','log')
set(gca, 'XTick', [10^0, 10^2, 10^4])
set(gca, 'YTick', [10^0, 10^2, 10^4])
set(gca, 'TickLength', [0.04, 0.01])
box on
hold off

nexttile
hold on
xlabel('\eta\it{F}\rm{^2_{cal}(b)}')
ylabel('\eta\it{F}\rm{^2_{obs}(b)}')
e1p5KHK0Het = errorbar(F2cal1p5KHK0Het, F2obs1p5KHK0Het, SigmaObs1p5KHK0Het, 'o', 'MarkerFaceColor', 'w');
e1p6K0KLHet = errorbar(F2cal1p6K0KLHet, F2obs1p6K0KLHet, SigmaObs1p6K0KLHet, 'o', 'MarkerFaceColor', 'w');
plot(LinPhHet, LinPhHet, 'k', 'LineWidth', 1)
legend([e1p5KHK0Het, e1p6K0KLHet], {'1.5 K, (HK0) Het.', '1.6 K, (0KL) Het.'});
axis([min(LinPhHet), max(LinPhHet), min(LinPhHet), max(LinPhHet)])
axis square
set(gca,'XScale','log','YScale','log')
set(gca, 'XTick', [10^0, 10^2, 10^4])
set(gca, 'YTick', [10^0, 10^2, 10^4])
set(gca, 'TickLength', [0.04, 0.01])
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Refinement\Refinement.eps';
exportgraphics(gcf, fdir, 'ContentType', 'vector')

%% Eu5In2Sb6 Susceptibility Analysis. Vincent Morano, 10/05/2022
% c-axis: 10/5/22 VSM. 1.59(2) mg, TBe343.
% b-axis: 10/24/22 MPMS. 17.17 mg, TBe339.
% a-axis: 04/15/23 MPMS. 16.06 mg, TBe339.

clear
close all

set(0, 'defaulttextinterpreter', 'tex')
set(0, 'defaultlegendinterpreter', 'tex')

transition7K = 7.2028; % Midpoint of sharp HC edge
transition14K = 14.0443; % Midpoint of sharp HC edge
massc = 1.59; % Sample mass in mg
massb = 17.17; % Sample mass in mg
massa = 16.06; % Sample mass in mg
molarMass = 5*(151.96)+2*(114.82)+6*(121.76); % Sample grams per mole formula unit
convFact1000c = 1/((massc/1e3)/molarMass)/5/1e3; % Convert from emu to emu/mol/Eu(/Oe).
convFact1000b = 1/((massb/1e3)/molarMass)/5/1e3; % Convert from emu to emu/mol/Eu(/Oe).
convFact1000a = 1/((massa/1e3)/molarMass)/5/1e3; % Convert from emu to emu/mol/Eu(/Oe).

% Import a-axis data
fileLoca = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MPMS\Eu5In2Sb6\TBe339\';
fileNamea = 'Eu5In2Sb6_TBe339_16p06mg_MPMS_MvT_aAxis_041523.dat';
headerLinesa = 45;
filea = strcat(fileLoca, fileNamea);
dataa = readtable(filea, 'FileType', 'delimitedtext', 'NumHeaderLines', headerLinesa, 'TrailingDelimitersRule', 'keep');
fida = fopen(filea); % Open file to take column headings as one line of text and edit
titlea = textscan(fida, '%s', 89, 'delimiter', ',', 'MultipleDelimsAsOne', 1, 'headerlines', headerLinesa-1);
fclose(fida);
titlea = titlea{1};
dataa.Properties.VariableNames = titlea; % Remove "Comment" from variable names

ind1000ZFCa = 2:448; % 1000 Oe, ZFC
ind1000FCa = 449:895; % 1000 Oe, FC

susc1000ZFCa = dataa(ind1000ZFCa, 'DC Moment Free Ctr (emu)').(1)*convFact1000a;
susc1000FCa = dataa(ind1000FCa, 'DC Moment Free Ctr (emu)').(1)*convFact1000a;

susc1000ZFCErra = dataa(ind1000ZFCa, 'DC Moment Err Free Ctr (emu)').(1)*convFact1000a;
susc1000FCErra = dataa(ind1000FCa, 'DC Moment Err Free Ctr (emu)').(1)*convFact1000a;

temp1000ZFCa = dataa(ind1000ZFCa, 'Temperature (K)').(1);
temp1000FCa = dataa(ind1000FCa, 'Temperature (K)').(1);

% Import b-axis data
fileLocb = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MPMS\Eu5In2Sb6\TBe339\';
fileNameb = 'Eu5In2Sb6_b_TBe339_17p17mg_b_102422.dat';
headerLinesb = 45;
fileb = strcat(fileLocb, fileNameb);
datab = readtable(fileb, 'FileType', 'delimitedtext', 'NumHeaderLines', headerLinesb, 'TrailingDelimitersRule', 'keep');
fidb = fopen(fileb); % Open file to take column headings as one line of text and edit
titleb = textscan(fidb, '%s', 89, 'delimiter', ',', 'MultipleDelimsAsOne', 1, 'headerlines', headerLinesb-1);
fclose(fidb);
titleb = titleb{1};
datab.Properties.VariableNames = titleb; % Remove "Comment" from variable names

ind1000ZFCb = 2:386; % 1000 Oe, ZFC
ind1000FCb = 388:1394; % 1000 Oe, FC

susc1000ZFCb = datab(ind1000ZFCb, 'DC Moment Free Ctr (emu)').(1)*convFact1000b;
susc1000FCb = datab(ind1000FCb, 'DC Moment Free Ctr (emu)').(1)*convFact1000b;

susc1000ZFCErrb = datab(ind1000ZFCb, 'DC Moment Err Free Ctr (emu)').(1)*convFact1000b;
susc1000FCErrb = datab(ind1000FCb, 'DC Moment Err Free Ctr (emu)').(1)*convFact1000b;

temp1000ZFCb = datab(ind1000ZFCb, 'Temperature (K)').(1);
temp1000FCb = datab(ind1000FCb, 'Temperature (K)').(1);

% Import c-axis data
fileLocc = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\VSM\Eu5In2Sb6\080522_TBe343_1p59mg_cAxis\';
fileNamec = 'MvT_1p59mg.dat';
headerLinesc = 23;
filec = strcat(fileLocc, fileNamec);
datac = readtable(filec, 'FileType', 'delimitedtext', 'NumHeaderLines', headerLinesc, 'TrailingDelimitersRule', 'keep');
fidc = fopen(filec); % Open file to take column headings as one line of text and edit
titlec = textscan(fidc, '%s', 57, 'delimiter', ',', 'MultipleDelimsAsOne', 1, 'headerlines', headerLinesc-1);
fclose(fidc);
titlec = titlec{1};
datac.Properties.VariableNames = titlec; % Remove "Comment" from variable names

ind1000ZFCc = 2:9724; % 1000 Oe, ZFC
ind1000FCc = 29129:38801; % 1000 Oe, FC. 9731:19288 had Issue where temperature was overshot and had to come back

susc1000ZFCc = datac(ind1000ZFCc, 'Moment (emu)').(1)*convFact1000c;
susc1000FCc = datac(ind1000FCc, 'Moment (emu)').(1)*convFact1000c;

susc1000ZFCErrc = datac(ind1000ZFCc, 'M. Std. Err. (emu)').(1)*convFact1000c;
susc1000FCErrc = datac(ind1000FCc, 'M. Std. Err. (emu)').(1)*convFact1000c;

temp1000ZFCc = datac(ind1000ZFCc, 'Temperature (K)').(1);
temp1000FCc = datac(ind1000FCc, 'Temperature (K)').(1);

% Plot susceptibility versus temperature
figure(1)
hold on
xlabel('\it{T}\rm{ (K)}', 'FontSize', 12)
ylabel('\chi (emu mol^{-1} Eu^{-1})', 'FontSize', 12)
% e1 = errorbar(temp1000ZFCa, susc1000ZFCa, susc1000ZFCErra, 'o', 'MarkerFaceColor', 'w');
% e2 = errorbar(temp1000FCa, susc1000FCa, susc1000FCErra, 'o', 'MarkerFaceColor', 'w');
% e3 = errorbar(temp1000ZFCb, susc1000ZFCb, susc1000ZFCErrb, 'o', 'MarkerFaceColor', 'w');
% e4 = errorbar(temp1000FCb, susc1000FCb, susc1000FCErrb, 'o', 'MarkerFaceColor', 'w');
% e5 = errorbar(temp1000ZFCc, susc1000ZFCc, susc1000ZFCErrc, 'o', 'MarkerFaceColor', 'w');
% e6 = errorbar(temp1000FCc, susc1000FCc, susc1000FCErrc, 'o', 'MarkerFaceColor', 'w');
e1 = plot(temp1000ZFCa, susc1000ZFCa, 'LineWidth', 2);
e2 = plot(temp1000FCa, susc1000FCa, 'LineWidth', 2);
e3 = plot(temp1000ZFCb, susc1000ZFCb, 'LineWidth', 2);
e4 = plot(temp1000FCb, susc1000FCb, 'LineWidth', 2);
e5 = plot(temp1000ZFCc, susc1000ZFCc, 'LineWidth', 2);
e6 = plot(temp1000FCc, susc1000FCc, 'LineWidth', 2);
%legend([e1, e2, e3, e4, e5, e6], {'a-axis ZFC', 'a-axis FC', 'b-axis ZFC', 'b-axis FC', 'c-axis ZFC', 'c-axis FC'}, 'fontsize', 10)
text(0.7, 0.9, '\it{H}\rm{ = 1000 Oe}', 'units', 'normalized', 'FontSize', 12)
text(0.32, 0.76, ['a (\color[rgb]{', num2str(e2.Color), '}FC\color{black}, \color[rgb]{', num2str(e1.Color), '}ZFC\color{black})'], 'units', 'normalized', 'FontSize', 11)
text(0.32, 0.39, ['b (\color[rgb]{', num2str(e4.Color), '}FC\color{black}, \color[rgb]{', num2str(e3.Color), '}ZFC\color{black})'], 'units', 'normalized', 'FontSize', 11)
text(0.17, 0.65, ['c (\color[rgb]{', num2str(e6.Color), '}FC\color{black}, \color[rgb]{', num2str(e5.Color), '}ZFC\color{black})'], 'units', 'normalized', 'FontSize', 11)
xlim([0, 35])
ylim([-0.05, 0.75])
set(gca, 'FontSize', 12)
pbaspect([16 9 1])
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\ChivT.eps';
%exportgraphics(gcf, fdir, 'ContentType', 'vector')

% Plot magnetization versus field
NA = 6.022e23;
% Import a-axis data
fileLocMvHa = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MPMS\Eu5In2Sb6\TBe339\';
fileNameMvHa = 'Eu5In2Sb6_TBe339_16p06mg_MPMS_MvH_10K_aAxis_041523.dat';
fileMvHa = strcat(fileLocMvHa, fileNameMvHa);
headerLinesMvHa = 45;
dataMvHa = readtable(fileMvHa, 'FileType', 'delimitedtext', 'NumHeaderLines', headerLinesMvHa, 'TrailingDelimitersRule', 'keep');
fidMvHa = fopen(fileMvHa); % Open file to take column headings as one line of text and edit
titleMvHa = textscan(fidMvHa, '%s', 89, 'delimiter', ',', 'MultipleDelimsAsOne', 1, 'headerlines', headerLinesMvHa-1);
fclose(fidMvHa);
titleMvHa = titleMvHa{1};
dataMvHa.Properties.VariableNames = titleMvHa; % Remove "Comment" from variable names
ind10KMvHa = 1:601; % 10 K, MvH, a
field10KMvHa = dataMvHa(ind10KMvHa, 'Magnetic Field (Oe)').(1);
convFact10KMvHa = 1/9.274e-21/(massa/1e3)*molarMass/NA/5; % Convert from emu to muB/Eu. See 03/30/2023 526 notes.
convFact10KMvHaemug = 1/(massa/1e3); % Convert from emu to emu/g.
convFact10KMvHaemumoleu = 1/((massa/1e3)/molarMass)/5; % Convert from emu to emu/mol/Eu.
mom10KMvHa = dataMvHa(ind10KMvHa, 'DC Moment Free Ctr (emu)').(1).*convFact10KMvHa;
mom10KMvHaemug = dataMvHa(ind10KMvHa, 'DC Moment Free Ctr (emu)').(1).*convFact10KMvHaemug;
mom10KMvHaemumoleu = dataMvHa(ind10KMvHa, 'DC Moment Free Ctr (emu)').(1).*convFact10KMvHaemumoleu;
mom10KMvHErra = dataMvHa(ind10KMvHa, 'DC Moment Err Free Ctr (emu)').(1).*convFact10KMvHa;
mom10KMvHErraemug = dataMvHa(ind10KMvHa, 'DC Moment Err Free Ctr (emu)').(1).*convFact10KMvHaemug;
mom10KMvHErraemumoleu = dataMvHa(ind10KMvHa, 'DC Moment Err Free Ctr (emu)').(1).*convFact10KMvHaemumoleu;
temp10KMvHa = dataMvHa(ind10KMvHa, 'Temperature (K)').(1);
figure(2)
hold on
xlabel('\it{H}\rm{ (Oe)}', 'FontSize', 12)
ylabel('\it{M}\rm{ (emu g^{-1})}', 'FontSize', 12)
e7 = errorbar(field10KMvHa, mom10KMvHaemug, mom10KMvHErraemug, 'o', 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 12)
pbaspect([16 9 1])
xlim([-2000, 2000])
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\MvH.eps';
%exportgraphics(gcf, fdir, 'ContentType', 'vector')

figure(3)
hold on
xlabel('\it{H}\rm{ (Oe)}', 'FontSize', 12)
ylabel('\it{M}\rm{ (\mu_B Eu^{-1})}', 'FontSize', 12)
e8 = errorbar(field10KMvHa, mom10KMvHa, mom10KMvHErra, 'o', 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 12)
pbaspect([16 9 1])
xlim([-2000, 2000])
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\MvHmuB.eps';
%exportgraphics(gcf, fdir, 'ContentType', 'vector')

figure(4)
hold on
xlabel('\it{H}\rm{ (Oe)}', 'FontSize', 12)
ylabel('\it{M}\rm{ (emu mol^{-1} Eu^{-1})}', 'FontSize', 12)
e7 = errorbar(field10KMvHa, mom10KMvHaemumoleu, mom10KMvHErraemumoleu, 'o', 'MarkerFaceColor', 'w');
set(gca, 'FontSize', 12)
pbaspect([16 9 1])
xlim([-2000, 2000])
ylim([-1100, 1100])
box on
hold off

fdir = 'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\Manuscript\MvHemumoleu.eps';
%exportgraphics(gcf, fdir, 'ContentType', 'vector')