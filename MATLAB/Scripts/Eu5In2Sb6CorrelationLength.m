% Vincent Morano, 9/8/2022
% Determine the correlation length for Eu5In2Sb6 SPINS and BT7 experiments

clear
close all

set(0, 'defaulttextinterpreter', 'latex')
set(0, 'defaultlegendinterpreter', 'latex')

% Fit to the (100) reflection from SPINS experiment with Voigt and pick out
% Lorentzian hwhm for correlation length
sqrt(8*log(2))
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
figure(1)
tiledlayout(2, 1, 'TileSpacing', 'compact')
nexttile
hold on
e1=errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
e2=errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
e3=errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%p1=plot(KCal(:,1), intCal(:,1), 'Color', 'r');
p2=plot(KCal(:,2), intCal(:,2), 'Color', 'b', 'LineWidth', 1);
p3=plot(KCal(:,3), intCal(:,3), 'Color', 'k', 'LineWidth', 1);
set(gca, 'FontSize', 12)
legend([e1, e2, e3],{'$T=20\ \mathrm{K}$', '$T=10\ \mathrm{K}$', '$T=1.5\ \mathrm{K}$'}, 'fontsize', 12)
text(0.2, 0.9, '$(100)$', 'Units', 'normalized', 'fontsize', 16)
text(0.2, 0.8, 'SPINS', 'Units', 'normalized', 'fontsize', 16)
%xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
xlim([min(K, [], 'all'), max(K, [], 'all')])
ylim([0 2200])
yticks(linspace(0, 1700, 5))
set(gca, 'TickLength',[0.02, 0.025]) % [2Dlength 3Dlength]
pbaspect([4 3 1])
box on
text(0.02, 0.9, '(a)', 'units', 'normalized', 'fontsize', 16)
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
        ylabel('$\chi^2_{\mathrm{r}}$')
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
    title([tempText{i}, ' $(100)$ Background: ', num2str(round(xFit(1, i), 4)), '$\pm$', num2str(round(xErr(1, i), 4))])
    xlabel('Background $\left(\frac{\mathrm{det.\ cts}}{\mathrm{7e5\ mon.\ cts}}\right)$')
    hold off
    subplot(2, 2, 2)
    hold on
    title([tempText{i}, ' $(100)$ Voigt Peak: ', num2str(round(xFit(2, i), 4)), '$\pm$', num2str(round(xErr(2, i), 4))])
    xlabel('Voigt Peak $\left(\frac{\mathrm{det.\ cts}}{\mathrm{7e5\ mon.\ cts}}\right)$')
    hold off
    subplot(2, 2, 3)
    hold on
    title([tempText{i}, ' $(100)$ Voigt Center: ', num2str(round(xFit(3, i), 4)), '$\pm$', num2str(round(xErr(3, i), 4))])
    xlabel('Voigt Center (deg.)')
    hold off
    subplot(2, 2, 4)
    hold on
    title([tempText{i}, ' $(100)$ Lorentzian FWHM: ', num2str(round(xFit(4, i), 5)), '$\pm$', num2str(round(xErr(4, i), 5))])
    xlabel('Lorentzian FWHM (deg.)')
    hold off
end

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\CorrBestSPINS.mat') % Save the Workspace

%% Add the (0,0,3/2) reflection from BT7 experiment
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
set(gca, 'FontSize', 12)
legend([e1, e2],{'$T=18\ \mathrm{K}$', '$T=1.6\ \mathrm{K}$'}, 'fontsize', 12)
text(0.2, 0.9, '$( 0 0 \frac{3}{2} )$', 'Units', 'normalized', 'fontsize', 16)
text(0.2, 0.8, 'BT7', 'Units', 'normalized', 'fontsize', 16)
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
xlim([min(K, [], 'all'), max(K, [], 'all')])
ylim([0 1600])
yticks(linspace(0, 1400, 5))
set(gca, 'TickLength',[0.02, 0.025]) % [2Dlength 3Dlength]
pbaspect([4 3 1])
box on
text(0.02, 0.9, '(b)', 'units', 'normalized', 'fontsize', 16)
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrBest.png')

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
    ylabel('$\chi^2_{\mathrm{r}}$')
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
title([tempText{i}, ' $\left( 0 0 \frac{3}{2} \right)$ Background: ', num2str(round(xFit(1, i), 4)), '$\pm$', num2str(round(xErr(1, i), 4))])
xlabel('Background $\left(\frac{\mathrm{det.\ cts}}{\mathrm{7e5\ mon.\ cts}}\right)$')
hold off
subplot(2, 2, 2)
hold on
title([tempText{i}, ' $\left( 0 0 \frac{3}{2} \right)$ Voigt Peak: ', num2str(round(xFit(2, i), 4)), '$\pm$', num2str(round(xErr(2, i), 4))])
xlabel('Voigt Peak $\left(\frac{\mathrm{det.\ cts}}{\mathrm{7e5\ mon.\ cts}}\right)$')
hold off
subplot(2, 2, 3)
hold on
title([tempText{i}, ' $\left( 0 0 \frac{3}{2} \right)$ Voigt Center: ', num2str(round(xFit(3, i), 4)), '$\pm$', num2str(round(xErr(3, i), 4))])
xlabel('Voigt Center (deg.)')
hold off
subplot(2, 2, 4)
hold on
title([tempText{i}, ' $\left( 0 0 \frac{3}{2} \right)$ Lorentzian FWHM: ', num2str(round(xFit(4, i), 5)), '$\pm$', num2str(round(xErr(4, i), 5))])
xlabel('Lorentzian FWHM (deg.)')
hold off

save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\CorrBestBT7.mat') % Save the Workspace

%% Test a few correlation lengths manually to see how shallow the gof min is: SPINS

clear
close all

set(0, 'defaulttextinterpreter', 'latex')
set(0, 'defaultlegendinterpreter', 'latex')

% Fit to the (100) reflection from SPINS experiment with Voigt and pick out
% Lorentzian hwhm for correlation length

% Variables for correlation length fit
gam = linspace(0.01, 0.2, 2e2); % Range of Lorentzian hwhm to plot
sigSPINS = 0.48742/(2*sqrt(2*log(2))); % See Neutron Analysis script
a = 12.51000; %A: angstroms ICSD_CollCode281238.cif
b = 14.58400; %B: angstroms
c = 4.62430; %C: angstroms
HKL = [1, 0, 0];
q = 2*pi*sqrt(HKL(1)^2/a^2 + HKL(2)^2/b^2 + HKL(3)^2/c^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
errPts = 0; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.3, 0.005]; % Determines range of iterated values for errorbar. [bg, peak, center]
offset = [0.0, 0.2, 0.005]; % Determines range of iterated values for errorbar
model = @(x, a3, sig, gam) x(1) + x(2).*voigt(a3, x(3), sig, gam); % Voigt function with flat background. center, x0, sigma, and gamma

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
dof = zeros(size(motor(1,:)));

% Fit for correlation length
for i = 1:length(motor(1,:))
    for j = 1:length(gam)
        tmp = sort(int(:,i)); % Sorted intensities
        bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
        cent0 = sum(motor(:,i).*(int(:,i)-bg0))./sum(int(:,i)-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
        peak0 = (max(int(:,i))-bg0)/voigt(cent0, cent0, sigSPINS, gam(j)); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
        x0 = [bg0, peak0, cent0];
        dof(i) = length(int(:,i))-length(x0); % Number of datapoints minus number of free parameters
        modelInput = @(x) model(x, motor(:,i), sigSPINS, gam(j));
        [xFit(:,i,j), redChi2Fit(:,i,j), xErr(:,i,j), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(int(:,i), intErr(:,i), modelInput, x0, errPts, fact, offset);
        
        fwhm(i,j) = 2*q*gam(j)*pi/180; % Fwhm from Lorentzian component of Voigt function
        corr(i,j) = 2/fwhm(i,j); % Factor of 2 is important to remember
    
        % Values for plotting the fit
        motorCal(:,i,j) = linspace(min(motor(:,i)), max(motor(:,i)), 500)';
        intCal(:,i,j) = model(xFit(:,i,j), motorCal(:,i), sigSPINS, gam(j));

        disp(['Finished loop ', num2str(length(gam)*(i-1)+j), ' out of ', num2str(length(motor(1,:))*length(gam)), '.' ])
    end

    % Convert to rlu, in HK0 looking at 100 with a3 scan so taking 0K0 projection
    qK(:,i) = q.*sind(mean(motor(:,i), 1).*ones(size(motor(:,i)))-motor(:,i)); % Here the a3 angle in real space doesn't correspond to the rotated angle in reciprocal space because the lattice is orthorhombic.
    K(:,i) = qK(:,i)./bStar;
    qKCal(:,i) = q.*sind(mean(motorCal(:,i), 1).*ones(size(motorCal(:,i)))-motorCal(:,i));
    KCal(:,i) = qKCal(:,i)./bStar;
end

% Plot the gof's versus correlation length
figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(2,1)
% nexttile
% hold on
% title('20 K')
% xlabel('$\xi (\mathrm{\AA})$')
% ylabel('$\chi^2_r$')
% yline(min(redChi2Fit(:,1,:))*(1+1/dof(1)), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
% scatter(corr(1,:), reshape(redChi2Fit(:,1,:), size(corr(1,:))), 'o')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
title('10 K')
xlabel('$\xi (\mathrm{\AA})$')
ylabel('$\chi^2_r$')
yline(min(redChi2Fit(:,2,:))*(1+1/dof(2)), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corr(2,:), reshape(redChi2Fit(:,2,:), size(corr(2,:))), 'o')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
title('1.5 K')
xlabel('$\xi (\mathrm{\AA})$')
ylabel('$\chi^2_r$')
yline(min(redChi2Fit(:,3,:))*(1+1/dof(3)), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corr(3,:), reshape(redChi2Fit(:,3,:), size(corr(3,:))), 'o')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrSPINS.png')

% Plot the best fits and the fits right next to the gof
corr20K = corr(1,:);
chiTrial20K = reshape(redChi2Fit(:,1,:), size(corr(1,:)));
chiUpper20K = min(chiTrial20K)*(1+1/dof(1));
[~, indBest20K] = min(chiTrial20K);
corrBest20K = corr20K(indBest20K);
minPt20K = min(abs(min(chiTrial20K)-chiTrial20K));
indMin20K = find(chiTrial20K<chiUpper20K & corr20K<corrBest20K, 1, 'last');
indMax20K = find(chiTrial20K<chiUpper20K & corr20K>corrBest20K, 1, 'first');

corr10K = corr(2,:);
chiTrial10K = reshape(redChi2Fit(:,2,:), size(corr(2,:)));
chiUpper10K = min(chiTrial10K)*(1+1/dof(2));
[~, indBest10K] = min(chiTrial10K);
corrBest10K = corr10K(indBest10K);
minPt10K = min(abs(min(chiTrial10K)-chiTrial10K));
indMin10K = find(chiTrial10K<chiUpper10K & corr10K<corrBest10K, 1, 'last');
indMax10K = find(chiTrial10K<chiUpper10K & corr10K>corrBest10K, 1, 'first');

corr1p5K = corr(3,:);
chiTrial1p5K = reshape(redChi2Fit(:,3,:), size(corr(3,:)));
chiUpper1p5K = min(chiTrial1p5K)*(1+1/dof(3));
[~, indBest1p5K] = min(chiTrial1p5K);
corrBest1p5K = corr1p5K(indBest1p5K);
minPt1p5K = min(abs(min(chiTrial1p5K)-chiTrial1p5K));
indMin1p5K = find(chiTrial1p5K<chiUpper1p5K & corr1p5K<corrBest1p5K, 1, 'last');
indMax1p5K = find(chiTrial1p5K<chiUpper1p5K & corr1p5K>corrBest1p5K, 1, 'first');

figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(2,2)
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% %plot(KCal(:,1), intCal(:,1,indMin20K), 'Color', 'b', 'LineWidth', 1);
% title('20 K Min')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% plot(KCal(:,1), intCal(:,1,indMax20K), 'Color', 'b', 'LineWidth', 1);
% title('20 K Max')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMin10K), 'Color', 'k', 'LineWidth', 1);
title('10 K Min')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMax10K), 'Color', 'k', 'LineWidth', 1);
title('10 K Max')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,3), intCal(:,3,indMin1p5K), 'Color', 'k', 'LineWidth', 1);
title('1.5 K Min')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,3), intCal(:,3,indMax1p5K), 'Color', 'k', 'LineWidth', 1);
title('1.5 K Max')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrSPINSPeaks.png')

% Plot the best fits and the fits right next to the gof using 10% gof
% threshold
chiUpper20K10Per = min(chiTrial20K)*1.1;
indMin20K10Per = find(chiTrial20K<chiUpper20K10Per & corr20K<corrBest20K, 1, 'last');
indMax20K10Per = find(chiTrial20K<chiUpper20K10Per & corr20K>corrBest20K, 1, 'first');

chiUpper10K10Per = min(chiTrial10K)*1.1;
indMin10K10Per = find(chiTrial10K<chiUpper10K10Per & corr10K<corrBest10K, 1, 'last');
indMax10K10Per = find(chiTrial10K<chiUpper10K10Per & corr10K>corrBest10K, 1, 'first');

chiUpper1p5K10Per = min(chiTrial1p5K)*1.1;
indMin1p5K10Per = find(chiTrial1p5K<chiUpper1p5K10Per & corr1p5K<corrBest1p5K, 1, 'last');
indMax1p5K10Per = find(chiTrial1p5K<chiUpper1p5K10Per & corr1p5K>corrBest1p5K, 1, 'first');

figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(2,2)
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% %plot(KCal(:,1), intCal(:,1,indMin20K10Per), 'Color', 'b', 'LineWidth', 1);
% title('20 K Min')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% plot(KCal(:,1), intCal(:,1,indMax20K10Per), 'Color', 'b', 'LineWidth', 1);
% title('20 K Max')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMin10K10Per), 'Color', 'k', 'LineWidth', 1);
title('10 K Min 10\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMax10K10Per), 'Color', 'k', 'LineWidth', 1);
title('10 K Max 10\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,3), intCal(:,3,indMin1p5K10Per), 'Color', 'k', 'LineWidth', 1);
title('1.5 K Min 10\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,3), intCal(:,3,indMax1p5K10Per), 'Color', 'k', 'LineWidth', 1);
title('1.5 K Max 10\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrSPINSPeaks10Per.png')

% Plot the best fits and the fits right next to the gof using 20% gof
% threshold
chiUpper20K20Per = min(chiTrial20K)*1.2;
indMin20K20Per = find(chiTrial20K<chiUpper20K20Per & corr20K<corrBest20K, 1, 'last');
indMax20K20Per = find(chiTrial20K<chiUpper20K20Per & corr20K>corrBest20K, 1, 'first');

chiUpper10K20Per = min(chiTrial10K)*1.2;
indMin10K20Per = find(chiTrial10K<chiUpper10K20Per & corr10K<corrBest10K, 1, 'last');
indMax10K20Per = find(chiTrial10K<chiUpper10K20Per & corr10K>corrBest10K, 1, 'first');

chiUpper1p5K20Per = min(chiTrial1p5K)*1.2;
indMin1p5K20Per = find(chiTrial1p5K<chiUpper1p5K20Per & corr1p5K<corrBest1p5K, 1, 'last');
indMax1p5K20Per = find(chiTrial1p5K<chiUpper1p5K20Per & corr1p5K>corrBest1p5K, 1, 'first');

figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(2,2)
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% %plot(KCal(:,1), intCal(:,1,indMin20K20Per), 'Color', 'b', 'LineWidth', 1);
% title('20 K Min')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
% plot(KCal(:,1), intCal(:,1,indMax20K20Per), 'Color', 'b', 'LineWidth', 1);
% title('20 K Max')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMin10K20Per), 'Color', 'k', 'LineWidth', 1);
title('10 K Min 20\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMax10K20Per), 'Color', 'k', 'LineWidth', 1);
title('10 K Max 20\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,3), intCal(:,3,indMin1p5K20Per), 'Color', 'k', 'LineWidth', 1);
title('1.5 K Min 20\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,3), int(:,3), intErr(:,3), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,3), intCal(:,3,indMax1p5K20Per), 'Color', 'k', 'LineWidth', 1);
title('1.5 K Max 20\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{7e5\ mon.\ cts} \right)$')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrSPINSPeaks20Per.png')
save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\CorrLengthSPINS.mat') % Save the Workspace

%% Test a few correlation lengths manually to see how shallow the gof min is: BT7

clear
close all

set(0, 'defaulttextinterpreter', 'latex')
set(0, 'defaultlegendinterpreter', 'latex')

% Fit to the (003/2) reflection from BT7 experiment with Voigt and pick out
% Lorentzian hwhm for correlation length

% Variables for correlation length fit
gam = linspace(0.001, 0.02, 2e2); % Range of Lorentzian hwhm to plot
sigBT7 = 0.38748/(2*sqrt(2*log(2))); % See Neutron Analysis script
a = 12.51000; %A: angstroms ICSD_CollCode281238.cif
b = 14.58400; %B: angstroms
c = 4.62430; %C: angstroms
HKL = [0, 0, 3/2];
q = 2*pi*sqrt(HKL(1)^2/a^2 + HKL(2)^2/b^2 + HKL(3)^2/c^2); % Define Q as ha*+kb*+lc*, calculated from definition of reciprocal lattice vectors, take square root of Q dotted with itself. Technically each datapoint has a different HKL, so averaging here.
errPts = 0; % Number of parameter iterations when calculating the errorbar
fact = [0.4, 0.3, 0.005]; % Determines range of iterated values for errorbar. [bg, peak, center]
offset = [0.0, 0.2, 0.005]; % Determines range of iterated values for errorbar
model = @(x, a3, sig, gam) x(1) + x(2).*voigt(a3, x(3), sig, gam); % Voigt function with flat background. center, x0, sigma, and gamma

% To convert to reciprocal lattice units later
aStar = 2.*pi./a;
bStar = 2.*pi./b;
cStar = 2.*pi./c;

% Import data
file1 = readcell(strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\fpx854193.bt7'), 'FileType', 'text', 'NumHeaderLines', 44, 'CommentStyle', '#', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');
file2 = readcell(strcat('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\NIST\Eu5In2Sb6\121120Data\fpx854147.bt7'), 'FileType', 'text', 'NumHeaderLines', 44, 'CommentStyle', '#', 'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join');

motor = [cell2mat(file1(:,13)), cell2mat(file2(:,13))];
int = [cell2mat(file1(:,10)), cell2mat(file2(:,10))];
intErr = sqrt(int);
dof = zeros(size(motor(1,:)));

% Fit for correlation length
for i = 1:length(motor(1,:))
    for j = 1:length(gam)
        tmp = sort(int(:,i)); % Sorted intensities
        bg0 = mean(tmp(1:3)); % Guess that the background is close to the average of the lowest couple of datapoints in the scan.
        cent0 = sum(motor(:,i).*(int(:,i)-bg0))./sum(int(:,i)-bg0); % Guess that the peak center is close to the average of the scanned angles weighted by their intensity.
        peak0 = (max(int(:,i))-bg0)/voigt(cent0, cent0, sigBT7, gam(j)); % Guess that the prefactor is the largest observed intensity in the scan minus the guessed background divided by the voigt function max since voigt(cent)!=1.
        x0 = [bg0, peak0, cent0];
        dof(i) = length(int(:,i))-length(x0); % Number of datapoints minus number of free parameters
        modelInput = @(x) model(x, motor(:,i), sigBT7, gam(j));
        [xFit(:,i,j), redChi2Fit(:,i,j), xErr(:,i,j), ~, ~, ~, ~, ~, ~, ~, ~] = fitRedChi2Err(int(:,i), intErr(:,i), modelInput, x0, errPts, fact, offset);
        
        fwhm(i,j) = 2*q*gam(j)*pi/180; % Fwhm from Lorentzian component of Voigt function
        corr(i,j) = 2/fwhm(i,j); % Factor of 2 is important to remember
    
        % Values for plotting the fit
        motorCal(:,i,j) = linspace(min(motor(:,i)), max(motor(:,i)), 500)';
        intCal(:,i,j) = model(xFit(:,i,j), motorCal(:,i), sigBT7, gam(j));

        disp(['Finished loop ', num2str(length(gam)*(i-1)+j), ' out of ', num2str(length(motor(1,:))*length(gam)), '.' ])
    end

    % Convert to rlu, in HK0 looking at 100 with a3 scan so taking 0K0 projection
    qK(:,i) = q.*sind(mean(motor(:,i), 1).*ones(size(motor(:,i)))-motor(:,i)); % Here the a3 angle in real space doesn't correspond to the rotated angle in reciprocal space because the lattice is orthorhombic.
    K(:,i) = qK(:,i)./bStar;
    qKCal(:,i) = q.*sind(mean(motorCal(:,i), 1).*ones(size(motorCal(:,i)))-motorCal(:,i));
    KCal(:,i) = qKCal(:,i)./bStar;
end

% Plot the gof's versus correlation length
figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
% tiledlayout(2,1)
% nexttile
% hold on
% title('18 K')
% xlabel('$\xi (\mathrm{\AA})$')
% ylabel('$\chi^2_r$')
% yline(min(redChi2Fit(:,1,:))*(1+1/dof(1)), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
% scatter(corr(1,:), reshape(redChi2Fit(:,1,:), size(corr(1,:))), 'o')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
hold on
title('1.6 K')
xlabel('$\xi (\mathrm{\AA})$')
ylabel('$\chi^2_r$')
yline(min(redChi2Fit(:,2,:))*(1+1/dof(2)), '--k', 'DisplayName', 'GOF Threshold', 'LineWidth', 1.0)
scatter(corr(2,:), reshape(redChi2Fit(:,2,:), size(corr(2,:))), 'o')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrBT7.png')

% Plot the best fits and the fits right next to the gof
corr18K = corr(1,:);
chiTrial18K = reshape(redChi2Fit(:,1,:), size(corr(1,:)));
chiUpper18K = min(chiTrial18K)*(1+1/dof(1));
[~, indBest18K] = min(chiTrial18K);
corrBest18K = corr18K(indBest18K);
minPt18K = min(abs(min(chiTrial18K)-chiTrial18K));
indMin18K = find(chiTrial18K<chiUpper18K & corr18K<corrBest18K, 1, 'last');
indMax18K = find(chiTrial18K<chiUpper18K & corr18K>corrBest18K, 1, 'first');

corr1p6K = corr(2,:);
chiTrial1p6K = reshape(redChi2Fit(:,2,:), size(corr(2,:)));
chiUpper1p6K = min(chiTrial1p6K)*(1+1/dof(2));
[~, indBest1p6K] = min(chiTrial1p6K);
corrBest1p6K = corr1p6K(indBest1p6K);
minPt1p6K = min(abs(min(chiTrial1p6K)-chiTrial1p6K));
indMin1p6K = find(chiTrial1p6K<chiUpper1p6K & corr1p6K<corrBest1p6K, 1, 'last');
indMax1p6K = find(chiTrial1p6K<chiUpper1p6K & corr1p6K>corrBest1p6K, 1, 'first');

figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(1,2)
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(KCal(:,1), intCal(:,1,indMin18K), 'Color', 'k', 'LineWidth', 1);
% title('18 K Min')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(KCal(:,1), intCal(:,1,indMax18K), 'Color', 'k', 'LineWidth', 1);
% title('18 K Max')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMin1p6K), 'Color', 'k', 'LineWidth', 1);
title('1.6 K Min')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMax1p6K), 'Color', 'k', 'LineWidth', 1);
title('1.6 K Max')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrBT7Peaks.png')

% Plot the best fits and the fits right next to the gof using 10% gof
% threshold
chiUpper18K10Per = min(chiTrial18K)*1.1;
indMin18K10Per = find(chiTrial18K<chiUpper18K10Per & corr18K<corrBest18K, 1, 'last');
indMax18K10Per = find(chiTrial18K<chiUpper18K10Per & corr18K>corrBest18K, 1, 'first');

chiUpper1p6K10Per = min(chiTrial1p6K)*1.1;
indMin1p6K10Per = find(chiTrial1p6K<chiUpper1p6K10Per & corr1p6K<corrBest1p6K, 1, 'last');
indMax1p6K10Per = find(chiTrial1p6K<chiUpper1p6K10Per & corr1p6K>corrBest1p6K, 1, 'first');

figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(1,2)
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(KCal(:,1), intCal(:,1,indMin18K), 'Color', 'k', 'LineWidth', 1);
% title('18 K Min')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(KCal(:,1), intCal(:,1,indMax18K), 'Color', 'k', 'LineWidth', 1);
% title('18 K Max')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMin1p6K10Per), 'Color', 'k', 'LineWidth', 1);
title('1.6 K Min 10\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMax1p6K10Per), 'Color', 'k', 'LineWidth', 1);
title('1.6 K Max 10\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrBT7Peaks10Per.png')

% Plot the best fits and the fits right next to the gof using 20% gof
% threshold
chiUpper18K20Per = min(chiTrial18K)*1.2;
indMin18K20Per = find(chiTrial18K<chiUpper18K20Per & corr18K<corrBest18K, 1, 'last');
indMax18K20Per = find(chiTrial18K<chiUpper18K20Per & corr18K>corrBest18K, 1, 'first');

chiUpper1p6K20Per = min(chiTrial1p6K)*1.2;
indMin1p6K20Per = find(chiTrial1p6K<chiUpper1p6K20Per & corr1p6K<corrBest1p6K, 1, 'last');
indMax1p6K20Per = find(chiTrial1p6K<chiUpper1p6K20Per & corr1p6K>corrBest1p6K, 1, 'first');

figure('Units', 'normalized', 'Position', [0.0, 0.3, 0.5, 0.6])
tiledlayout(1,2)
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(KCal(:,1), intCal(:,1,indMin18K), 'Color', 'k', 'LineWidth', 1);
% title('18 K Min')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
% pbaspect([16 9 1])
% box on
% hold off
% 
% nexttile
% hold on
% errorbar(K(:,1), int(:,1), intErr(:,1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(KCal(:,1), intCal(:,1,indMax18K), 'Color', 'k', 'LineWidth', 1);
% title('18 K Max')
% xlabel('(0K0) (r.l.u.)')
% ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
% pbaspect([16 9 1])
% box on
% hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMin1p6K20Per), 'Color', 'k', 'LineWidth', 1);
title('1.6 K Min 20\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
pbaspect([16 9 1])
box on
hold off

nexttile
hold on
errorbar(K(:,2), int(:,2), intErr(:,2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(KCal(:,2), intCal(:,2,indMax1p6K20Per), 'Color', 'k', 'LineWidth', 1);
title('1.6 K Max 20\%')
xlabel('(0K0) (r.l.u.)')
ylabel('Intensity $\left( \mathrm{det.\ cts} / \mathrm{15\ sec.} \right)$')
pbaspect([16 9 1])
box on
hold off

saveas(gcf,'C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\Figures\Eu5In2Sb6\CorrelationLength\CorrBT7Peaks20Per.png')
save('C:\Users\Vincent Morano\OneDrive - Johns Hopkins University\Lab\MATLAB\MAT\Eu5In2Sb6\CorrLengthBT7.mat') % Save the Workspace