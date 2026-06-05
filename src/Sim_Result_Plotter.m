% Plotting and Analyzing Simulation Results for All Lesion Sizes
clear; clc; close all;
load("results_keyhole.mat")

% --- Define the lesion sizes to iterate through ---
lesion_fields = {'Measured_Contrast_L', 'Measured_Contrast_M', 'Measured_Contrast_S', 'Measured_Contrast_XS'};
lesion_names = {'2 cm', '1 cm', '5 mm', '2.5 mm'};
field_suffixes = {'L', 'M', 'S', 'XS'}; % For valid dynamic struct field naming (no spaces)

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
        
        % Write back to the struct dynamically 
        resultsStruct2_flat(ii).(['Time_of_Peak_' field_suffixes{s}]) = t_peak;
        resultsStruct2_flat(ii).(['Peak_Enhancement_' field_suffixes{s}]) = c_max;
        
        % Calculate Euclidean Deviation from ground truth peak
        all_dev(ii, s) = sqrt((t_peak - gt_tp)^2 + (c_max - gt_max)^2);
    end
end

% --- 2.5 Sort the struct and arrays natively ---
[~, sort_idx] = sort([resultsStruct2_flat.Time_of_Peak_XS], 'ascend');
resultsStruct2_flat = resultsStruct2_flat(sort_idx);

% Apply the exact same sort to our pre-allocated arrays
all_tp = all_tp(sort_idx, :);
all_cmax = all_cmax(sort_idx, :);
all_dev = all_dev(sort_idx, :);

% Safely extract parameters for plotting
all_pA = [resultsStruct2_flat.pA]';
all_pB = [resultsStruct2_flat.pB]'; 


%% --- 3. Plot Peak Enhancement vs Time-to-Peak ---
fig2 = figure('Name', 'Enhancement vs Time-to-Peak', 'Position', [150, 150, 1200, 800]);
t_peak = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t_peak, 'Peak Enhancement vs Time-to-Peak by Lesion Size', 'FontWeight', 'bold', 'FontSize', 14);

% =========================================================================
% TOGGLE AXIS BEHAVIOR HERE:
% true  = All subplots share the exact same X and Y limits
% false = Each subplot scales independently to its own local maxima
use_global_axes = true; 
% =========================================================================

% Define Global Axis Limits (calculated once just in case)
global_xlim = [0, 2.5]; 
global_ylim = [0, 2.5];

ax_array = gobjects(1, 4); 

% Define distinct marker shapes for the different pB values
marker_shapes = {'o', 'v', '^', 'd', 'v', 'p', 'h'};
unique_pB = unique(all_pB);



for s = 1:4
    ax = nexttile(t_peak);
    ax_array(s) = ax; 
    hold(ax, 'on');
    
    % Loop through each unique pB value and plot its subset of data
    for b_idx = 1:length(unique_pB)
        current_pB = unique_pB(b_idx);
        idx = (all_pB == current_pB); % Find which rows match this pB
        
        m_shape = marker_shapes{mod(b_idx - 1, length(marker_shapes)) + 1};
        
        scatter(ax, all_tp(idx, s), all_cmax(idx, s), 45, all_pA(idx), ...
            m_shape, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
            'DisplayName', sprintf('pB = %.2f', current_pB));
    end
    
    % Reference lines for ground truth
    yline(ax, gt_max, 'r-', 'LineWidth', 1.5, 'DisplayName', 'GT Max Enhancement');
    %xline(ax, gt_tp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'GT Time-to-Peak');
    
    grid(ax, 'on');
    title(ax, lesion_names{s});
    xlabel(ax, 'Time of Peak Enhancement (s)');
    ylabel(ax, 'Peak Enhancement Value');
    
    % --- Apply limits based on user toggle ---
    if use_global_axes
        xlim(ax, global_xlim);
        ylim(ax, global_ylim);
    else
        % Calculate limits dynamically for THIS subplot only
        local_xlim = [0, max([all_tp(:, s); gt_tp]) * 1.1];
        local_ylim = [0, max([all_cmax(:, s); gt_max]) * 1.1];
        xlim(ax, local_xlim);
        ylim(ax, local_ylim);
    end
    
    if s == 2
        legend(ax, 'Location', 'bestoutside');
    end
    hold(ax, 'off');
end

% --- Link axes and add Colorbar ---
if use_global_axes
    linkaxes(ax_array, 'xy');
else
    linkaxes(ax_array, 'off'); % Ensures they act independently if local axes are used
end

colormap(fig2, 'parula'); 
cb = colorbar(ax); 
cb.Layout.Tile = 'east'; 
cb.Label.String = 'pA Value';
cb.Label.FontWeight = 'bold';


%% --- 4. Plot Deviation from Ground Truth ---
% We sort the x-axis index by the Time-to-Peak of the LARGE lesion 
% so all 4 subplots share a consistent, meaningful X-axis progression.
[~, sort_idx] = sort([resultsStruct2_flat.pA], 'ascend');
all_pA = sort(all_pA,'ascend');
[~, sort_idx_tp] = sort(all_tp(:, 1), 'ascend');

fig3 = figure('Name', 'Deviation from Ground Truth', 'Position', [200, 200, 1200, 800]);
t_dev = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(t_dev, 'Euclidean Distance from Ground Truth Kinetics', 'FontWeight', 'bold', 'FontSize', 14);

for s = 1:4
    ax = nexttile(t_dev);
    hold(ax, 'on');

    plot(ax, all_pA, all_dev(sort_idx_tp, s), '-o', 'LineWidth', 1.5, 'Color', '#D95319');

    grid(ax, 'on');
    title(ax, lesion_names{s});
    xlabel(ax, 'pA');
    ylabel(ax, 'Euclidean Distance');
    yticklabels([]);
    hold(ax, 'off');
end