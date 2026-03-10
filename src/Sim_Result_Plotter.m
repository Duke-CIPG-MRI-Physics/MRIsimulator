
%Plotting Simulation Results
clear;clc;close all;
load("Simulation Results.mat")

results = Simulation_Results.Results;

for ii = 1:length(results)
    tp_ii = results(ii).Timepoints;
    Simulation_Results.Results(ii).Temporal_Resolution = (tp_ii(end)-tp_ii(end-1));
end

results = struct2table(Simulation_Results.Results);
results = sortrows(results,"Temporal_Resolution","ascend");
results = table2struct(results);


figure;
t = tiledlayout(1, 1, 'TileSpacing', 'compact');
ax = nexttile(t);
hold(ax, 'on');

% Plot the Single Simulated Ground Truth Curve FIRST
plot(ax, results(1).Sim_Timepoints, results(1).Sim_Contrast,...
    'LineWidth', 3, 'DisplayName', 'Analytical Ground Truth');

% Define a color order for the measured curves
colors = lines(11);

% Loop through and overlay the measured TWIST data
for ii = 1:11
    meas_name = sprintf('TWIST (pA=%.2f, pB=%.2f)', results(ii).pA, results(ii).pB);
    
    plot(ax, results(ii).Timepoints, results(ii).Measured_Contrast, ...
        '-o', 'Color', colors(ii,:), 'LineWidth', .5, 'DisplayName', meas_name);
end

hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'Time Since Injection (s)', 'FontWeight', 'bold');
ylabel(ax, 'ROI Contrast Intensity', 'FontWeight', 'bold');
title(ax, 'TWIST Kinetics: Reconstructed vs. True Enhancement');
legend(ax, 'Location', 'bestoutside');