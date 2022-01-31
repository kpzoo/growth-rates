% Compare growth rates and reproduction numbers
clearvars; clc; close all; tic;

% Assumptions and notes
% - estimate Rt with EpiFilter and rt with Wallinga/non parametric

% Directory and where saving
thisDir = cd; saveFol = 'Results/'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Setup parameters and simulate epidemics 

% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;
% Choose scenarios and fix serial interval
scenNo = 6; distNo = 2;

% Simulation parameters and warning
Iwarn = 1; simVals = setupScenario(scenNo);
simVals.offset = 0;
% Simulate epidemic scenarios and truncate
while Iwarn
    [Iday, Lday, Rtrue, tday, Iwarn, distvals] = epiSimScenOffset(scenNo,...
        nday0, distNo, simVals);
end
if Iwarn
    warning('Sequences of zero incidence');
end

% Truncated observation period 
nday = length(tday);
% Saving data and figs
namstr = [num2str(nday) '_' num2str(scenNo) '_' num2str(distNo) '_'];

% Log derivative with removal of delays
rlam = diff(log(Lday)); rdelay = round(distvals.omega/2);

% Gamma distribution shape and scale
shape = distvals.pm; scale = distvals.omega/shape;
% True growth rate
rtrue = (Rtrue.^(1/shape) - 1)/scale;

%% Estimate epidemic transmissibility

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% EpiFilter estimate and prediction 
[Re, Ip, ~, re] = allFilSmoothGrow(Rgrid, m, eta, nday, p0, Lday, Iday, distvals);

% Other estimates by smoothing incidence
Iemp = smoothdata(Iday, 'sgolay', 50, 'Degree', 3); 
% Smooth and adjust for delay
remp = diff(log(Iemp)); remp = smoothdata(remp, 'sgolay', 90, 'Degree', 3);

% Another estimate from mean predicted incidence
Ipred = Ip.mean(:,2); Ipred = smoothdata(Ipred, 'sgolay', 50, 'Degree', 3); 
rpred = diff(log(Ipred)); rpred = smoothdata(rpred, 'sgolay', 40, 'Degree', 3);

%% Visualise results

% Estimate of Rt
figure;
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tday', Re.mean(:, 2), Re.low(:, 2), Re.high(:, 2), 'r');
plot(tday(2:end), ones(1, nday-1), '--', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('$\hat{R}_{t}$', 'FontSize', 20);
xlabel('time, $t$ (days)', 'FontSize', 20);
xlim([tday(2) tday(end)]); 

% Incidence prediction
figure;
scatter(tday(2:end)', Iday(2:end), 'Marker', 'o', 'SizeData', 50, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
hold on;
plotCIRaw(tday(2:end)', Ip.mean(:, 2), Ip.low(:, 2), Ip.high(:, 2), 'b');
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('$\hat{I}_{t}$', 'FontSize', 20);
xlabel('time, $t$ (days)', 'FontSize', 20);
xlim([tday(2) tday(end)]);

% Compare naive growth rates with R-based ones 
figure;
plot(tday(1:end-rdelay)', rlam(rdelay:end), 'Color', 'b', 'LineWidth', 2);
hold on;
plot(tday(2+rdelay:end)', remp(1:end-rdelay), 'Color', grey1, 'LineWidth', 2);
plotCIRaw(tday', re.mean(:, 2), re.low(:, 2), re.high(:, 2), 'r');
plot(tday(2:end), zeros(1, nday-1), '--', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('$\hat{r}_{t}$', 'FontSize', 20);
xlabel('time, $t$ (days)', 'FontSize', 20);
xlim([tday(2) tday(end)]);

% Single figure with both Rt and rt
figure;
subplot(2, 1, 1);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tday', Re.mean(:, 2), Re.low(:, 2), Re.high(:, 2), 'r');
plot(tday(2:end), ones(1, nday-1), '--', 'Color', grey1, 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('$\hat{R}_{t}$', 'FontSize', 20);
xlabel('time, $t$ (days)', 'FontSize', 20);
xlim([tday(2) tday(end)]);
subplot(2, 1, 2);
plot(tday, rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
plot(tday(1:end-rdelay)', rlam(rdelay:end), 'Color', 'b', 'LineWidth', 2);
plotCIRaw(tday', re.mean(:, 2), re.low(:, 2), re.high(:, 2), 'r');
plot(tday(2:end), zeros(1, nday-1), '--', 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('$\hat{r}_{t}$', 'FontSize', 20);
xlabel('time, $t$ (days)', 'FontSize', 20);

%% Publishable

% Main figure with rt, Rt, It, It smoothed
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
subplot(2, 2, 1);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
h = plotCIRaw(tday', Re.mean(:, 2), Re.low(:, 2), Re.high(:, 2), 'r');
plot(tday(2:end), ones(1, nday-1), '--', 'Color', 'k', 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('reproduction number', 'FontSize', 20);
%leg = legend('$R_t$', '', '$\hat{R}_t$', 'Location', 'Best');
%set(leg,'Box','off')
title('A');

subplot(2, 2, 3);
plotCIRaw(tday', re.mean(:, 2), re.low(:, 2), re.high(:, 2), 'r');
hold on;
plot(tday(1:end-rdelay)', rlam(rdelay:end), 'Color', 'b', 'LineWidth', 2);
plot(tday(2+rdelay:end)', remp(1:end-rdelay), 'Color', grey2, 'LineWidth', 2);
plot(tday(2:end), zeros(1, nday-1), '--', 'Color', 'k', 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('growth rate', 'FontSize', 20);
xlabel('time, $t$', 'FontSize', 20);
title('C');

subplot(2, 2, 2);
scatter(tday(2:end)', Iday(2:end), 'Marker', 'o', 'SizeData', 40, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.6);
%stairs(tday(2:end)', Ideme(2:end), 'k', 'LineWidth', 2);
hold on;
plotCIRaw(tday(2:end)', Ip.mean(:, 2), Ip.low(:, 2), Ip.high(:, 2), 'r');
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('case incidence', 'FontSize', 20);
title('B');

subplot(2, 2, 4);
%stairs(tday(2:end)', Ideme(2:end), 'Color', grey1, 'LineWidth', 2);
scatter(tday(2:end)', Iday(2:end), 'Marker', 'o', 'SizeData', 40, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.6);
hold on;
plot(tday(1:end-2*rdelay+1), Lday(2*rdelay:end), 'b', 'LineWidth', 2);
plot(tday, Iemp, 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('smoothed incidence', 'FontSize', 20);
xlabel('time, $t$', 'FontSize', 20);
title('D');

if saveFig
    cd(saveFol);
    saveas(gcf, ['Fig1a' namstr 'grow'], 'fig');
    cd(thisDir);
end


%% Publishable with true rt

% Main figure with rt, Rt, It, It smoothed
figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
subplot(2, 2, 1);
plot(tday, Rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
h = plotCIRaw(tday', Re.mean(:, 2), Re.low(:, 2), Re.high(:, 2), 'r');
plot(tday(2:end), ones(1, nday-1), '--', 'Color', 'k', 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('reproduction number', 'FontSize', 20);
%leg = legend('$R_t$', '', '$\hat{R}_t$', 'Location', 'Best');
%set(leg,'Box','off')
title('A');

h(1) = subplot(2, 2, 3);
plot(tday', re.mean(:, 2), 'r', 'LineWidth', 2);
hold on;
plot(tday(1:end-rdelay)', rlam(rdelay:end), 'Color', 'b', 'LineWidth', 2);
plot(tday(2+rdelay:end)', remp(1:end-rdelay), 'Color', grey2, 'LineWidth', 2);
plot(tday(2:end), zeros(1, nday-1), '--', 'Color', 'k', 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('growth rate', 'FontSize', 20);
xlabel('time, $t$', 'FontSize', 20);
title('C');

h(2) = subplot(2, 2, 2);
plot(tday, rtrue, '-', 'color', 'k', 'linewidth', 2);
hold on;
plotCIRaw(tday', re.mean(:, 2), re.low(:, 2), re.high(:, 2), 'r');
%plot(tday(1:end-rdelay)', rlam(rdelay:end), 'Color', 'b', 'LineWidth', 2);
%plot(tday(2+rdelay:end)', remp(1:end-rdelay), 'Color', grey2, 'LineWidth', 2);
plot(tday(2:end), zeros(1, nday-1), '--', 'Color', 'k', 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('growth rate', 'FontSize', 20);
xlabel('time, $t$', 'FontSize', 20);
title('B'); linkaxes(h, 'xy');


subplot(2, 2, 4);
%stairs(tday(2:end)', Ideme(2:end), 'Color', grey1, 'LineWidth', 2);
scatter(tday(2:end)', Iday(2:end), 'Marker', 'o', 'SizeData', 40, 'MarkerFaceColor',...
    grey1, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.6);
hold on;
plot(tday(1:end-2*rdelay+1), Lday(2*rdelay:end), 'b', 'LineWidth', 2);
plot(tday, Iemp, 'Color', grey2, 'LineWidth', 2);
grid off; box off; hold off; xlim([tday(2) tday(end)]);
ylabel('smoothed incidence', 'FontSize', 20);
xlabel('time, $t$', 'FontSize', 20);
title('D');

if saveFig
    cd(saveFol);
    saveas(gcf, ['Fig1b' namstr 'grow'], 'fig');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Remove unneeded variables and save
if saveTrue
    cd(saveFol);
    save(['growth' namstr  'data' '.mat']);
    cd(thisDir);
end
