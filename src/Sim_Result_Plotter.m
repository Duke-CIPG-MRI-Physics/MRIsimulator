% Plotting and Analyzing Simulation Results for All Lesion Sizes
clear; clc; close all;
load("Results_test.mat")

% --- Define the lesion sizes to iterate through ---
lesion_fields = {'Measured_Contrast_L', 'Measured_Contrast_M', 'Measured_Contrast_S', 'Measured_Contrast_XS'};
lesion_names = {'Large Lesion', 'Medium Lesion', 'Small Lesion', 'Extra Small Lesion'};
size_labels = {'L', 'M', 'S', 'XS'}; % For dynamic struct field naming

% Find actual Ground Truth metrics once
[gt_max, gt_idx] = max(Simulated.contrast);
gt_tp = Simulated.timepoints(gt_idx);

num_results = length(resultsStruct2_flat);
num_to_plot = min(10, num_results); % Plot up to 10 lines for visual clarity
colors = lines(num_to_plot);

% Pre-allocate arrays for metrics across all simulations and sizes (Sims x Sizes)
all_tp = zeros(num_results, 4);
all_cmax = zeros(num_results, 4);
all_dev = zeros(num_results, 4);

%% --- 1. Plot Kinetics for All 4 Sizes ---
fig1 = figure('Name', 'TWIST Kinetics', 'Position', [100, 100, 1200, 800]);
t_kinetics = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t_kinetics, 'TWIST Kinetics: Reconstructed vs. True Enhancement', 'FontWeight', 'bold', 'FontSize', 14);

for s = 1:4
    ax = nexttile(t_kinetics);
    hold(ax, 'on');
    
    % Plot Ground Truth
    plot(ax, Simulated.timepoints, Simulated.contrast, 'k--', 'LineWidth', 2, 'DisplayName', 'Analytical GT');
    
    % Loop through and overlay the measured TWIST data
    for ii = 1:num_to_plot
        meas_name = sprintf('pA=%.2f, pB=%.2f', resultsStruct2_flat(ii).pA, resultsStruct2_flat(ii).pB);
        plot(ax, resultsStruct2_flat(ii).Timepoints, resultsStruct2_flat(ii).(lesion_fields{s}), ...
            '-o', 'Color', colors(ii,:), 'LineWidth', 0.5, 'MarkerSize', 3, 'DisplayName', meas_name);
    end
    
    grid(ax, 'on');
    title(ax, lesion_names{s});
    xlabel(ax, 'Time Since Injection (s)');
    ylabel(ax, 'ROI Contrast Intensity');
    if s == 2 % Only put the legend on the top right tile to save space
        legend(ax, 'Location', 'bestoutside');
    end
    hold(ax, 'off');
end

%% --- 2. Calculate Metrics Across All Sizes ---
for ii = 1:num_results
    for s = 1:4
        % Extract Max and Time-to-Peak for the specific size
        [c_max, tp_idx] = max(resultsStruct2_flat(ii).(lesion_fields{s}));
        t_peak = resultsStruct2_flat(ii).Timepoints(tp_idx);
        
        % Store in arrays for plotting later
        all_tp(ii, s) = t_peak;
        all_cmax(ii, s) = c_max;
        
        % Write back to the struct dynamically (e.g., .Time_of_Peak_L)
        resultsStruct2_flat(ii).(['Time_of_Peak_' size_labels{s}]) = t_peak;
        resultsStruct2_flat(ii).(['Peak_Enhancement_' size_labels{s}]) = c_max;
        
        % Calculate Euclidean Deviation from ground truth peak
        all_dev(ii, s) = sqrt((t_peak - gt_tp)^2 + (c_max - gt_max)^2);
    end
end

% --- 2.5 Sort the struct natively (Bypasses struct2table errors) ---
% Extract all the Time_of_Peak_L values into an array, sort them, and use the index
[~, sort_idx] = sort([resultsStruct2_flat.Time_of_Peak_L], 'ascend');
resultsStruct2_flat = resultsStruct2_flat(sort_idx);

%% --- 3. Plot Peak Enhancement vs Time-to-Peak ---
fig2 = figure('Name', 'Enhancement vs Time-to-Peak', 'Position', [150, 150, 1200, 800]);
t_peak = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t_peak, 'Peak Enhancement vs Time-to-Peak by Lesion Size', 'FontWeight', 'bold', 'FontSize', 14);

for s = 1:4
    ax = nexttile(t_peak);
    hold(ax, 'on');
    
    scatter(ax, all_tp(:, s), all_cmax(:, s), 30, 'filled', 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Simulations');
    
    % Reference lines for ground truth
    yline(ax, gt_max, 'k--', 'LineWidth', 1.5, 'DisplayName', 'GT Max Enhancement');
    xline(ax, gt_tp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'GT Time-to-Peak');
    
    grid(ax, 'on');
    title(ax, lesion_names{s});
    xlabel(ax, 'Time of Peak Enhancement (s)');
    ylabel(ax, 'Peak Enhancement Value');
    if s == 2
        legend(ax, 'Location', 'bestoutside');
    end
    hold(ax, 'off');
end

%% --- 4. Plot Deviation from Ground Truth ---
% We sort the x-axis index by the Time-to-Peak of the LARGE lesion 
% so all 4 subplots share a consistent, meaningful X-axis progression.
[~, sort_idx] = sort(all_tp(:, 1), 'ascend');

fig3 = figure('Name', 'Deviation from Ground Truth', 'Position', [200, 200, 1200, 800]);
t_dev = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t_dev, 'Euclidean Distance from Ground Truth Kinetics', 'FontWeight', 'bold', 'FontSize', 14);

for s = 1:4
    ax = nexttile(t_dev);
    hold(ax, 'on');
    
    plot(ax, all_dev(sort_idx, s), '-o', 'LineWidth', 1.5, 'Color', '#D95319');
    
    grid(ax, 'on');
    title(ax, lesion_names{s});
    xlabel(ax, 'Simulation Index (Sorted by Large Lesion Wash-in)');
    ylabel(ax, 'Euclidean Distance');
    hold(ax, 'off');
end