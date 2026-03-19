
%Plotting Simulation Results
clear;clc;close all;
load("Results_test.mat")
load("ground_truth_contrast.mat")

% for ii = 1:length(results)
%     tp_ii = results(ii).Timepoints;
%     Simulation_Results.Results(ii).Temporal_Resolution = (tp_ii(end)-tp_ii(end-1));
% end
% 
% results = struct2table(Simulation_Results.Results);
% results = sortrows(results,"Temporal_Resolution","ascend");
% results = table2struct(results);


figure;
t = tiledlayout(1, 1, 'TileSpacing', 'compact');
ax = nexttile(t);
hold(ax, 'on');

% Plot the Single Simulated Ground Truth Curve FIRST
plot(ax, sim_tp, sim_contrast,...
    'LineWidth', 3, 'DisplayName', 'Analytical Ground Truth');

% Define a color order for the measured curves
colors = lines(10);
colors(11,:) = [1 0 0];

% Loop through and overlay the measured TWIST data
for ii = 1:10:20
    meas_name = sprintf('TWIST (pA=%.2f, pB=%.2f)', resultsStruct2_flat(ii).pA, resultsStruct2_flat(ii).pB);
    
    plot(ax, resultsStruct2_flat(ii).Timepoints, resultsStruct2_flat(ii).Measured_Contrast, ...
        '-o',  'LineWidth', .5, 'DisplayName', meas_name);
end

hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'Time Since Injection (s)', 'FontWeight', 'bold');
ylabel(ax, 'ROI Contrast Intensity', 'FontWeight', 'bold');
title(ax, 'TWIST Kinetics: Reconstructed vs. True Enhancement');
legend(ax, 'Location', 'bestoutside');

%% - Calculate time of peak enhancement for each sim settings

for ii = 1:length(resultsStruct2_flat)
   [contrast_max(ii),tp_idx(ii)]  = max(resultsStruct2_flat(ii).Measured_Contrast);
   tp(ii) = resultsStruct2_flat(ii).Timepoints(tp_idx(ii));

end


fieldName = "Time_of_Peak_Enhancement";

for i = 1:numel(resultsStruct2_flat)
    resultsStruct2_flat(i).(fieldName) = tp(i);
end

fieldName = "Peak_Enhancement";

for i = 1:numel(resultsStruct2_flat)
    resultsStruct2_flat(i).(fieldName) = contrast_max(i);
end


results = struct2table(resultsStruct2_flat);
results = sortrows(results,"Time_of_Peak_Enhancement","ascend");
results = table2struct(results);

figure
scatter(tp,contrast_max)
xlabel("Time of Peak Enhancement, Relative to Ground Truth");
ylabel("Peak Enhancement Value")
title("Enhancement Value vs Speed")

hold on
t = 0:.1:45;
v = ones(size(t)) *max(sim_contrast);
plot(t,v,"LineWidth",3);
legend("Simulations","Ground Truth Enhancement Value")
hold off

%% -- Calculate deviation from ground truth peak enhancement

deviation_from_gt_peak = sqrt((tp-0).^2 + (contrast_max-max(sim_contrast)).^2);
figure
plot(deviation_from_gt_peak)