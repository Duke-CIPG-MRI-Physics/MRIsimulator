clear;clc;close all;

% Inputs
Matrix_Size_Acquired = [1, 100, 100]; % [Freq, Phase, Slice] %frequency doesn't matter for this demo
FOV_acquired = [1, 100, 300];
pA = 0.05;
pB = .1; 
N_Measurements = 10; 
TR = 5e-3;
R = [1, 1];
PF_Factor = [1, 1]; 
orderingOptions = getTWISTOrderingOptions(struct( ...
    'radialBinWidthMode', "max", ...
    'bSubsetAssignment', "contiguous"));

[Complete_Sampling_Table, TWIST_Stats] = ...
    Ultrafast_Sampling( ...
        Matrix_Size_Acquired, FOV_acquired, pA, pB, N_Measurements, TR, R, PF_Factor, ...
        "forward", "single_anchor", orderingOptions);

[regionA, phaseEncodeTable] = getRegionA( ...
    Matrix_Size_Acquired, FOV_acquired, pA, PF_Factor, R, orderingOptions);

sz_4D = [Matrix_Size_Acquired, N_Measurements + 1];

sampled_mask = false(sz_4D);

sampled_mask(Complete_Sampling_Table.("Linear Index")) = true;


visualization_matrix = squeeze(sampled_mask(1, :, :, :));

imshow3D(visualization_matrix);

radialBinGrid = zeros(Matrix_Size_Acquired(2), Matrix_Size_Acquired(3));
radialBinGrid(phaseEncodeTable.LinearIndex) = phaseEncodeTable.RadialBin;

figure('Name', 'Ultrafast Sampling Diagnostics', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(regionA');
axis image;
set(gca, 'YDir', 'normal');
title('Region A');
xlabel('Phase Encode');
ylabel('Slice Encode');

nexttile;
imagesc(radialBinGrid');
axis image;
set(gca, 'YDir', 'normal');
title(sprintf('Radial Bins (%s dk)', orderingOptions.radialBinWidthMode));
xlabel('Phase Encode');
ylabel('Slice Encode');
colorbar;

plotTWISTFrameOrdering( ...
    Complete_Sampling_Table, Matrix_Size_Acquired, regionA, ...
    "Ultrafast B-Frame Ordering", 9);

%% --- Multi-Color Trajectory Animation (Fixed Speed)
% 1. Filter: Use only unique Phase-Slice events
traj_table = Complete_Sampling_Table(Complete_Sampling_Table.Frequency == 1, :);

% 2. Setup Figure
fig = figure('Name', 'Multi-Color Trajectory Animation', 'Color', 'k'); 
ax = axes('Parent', fig);
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w');
hold(ax, 'on');
axis(ax, 'square');
grid(ax, 'on');

xlim(ax, [0, Matrix_Size_Acquired(2) + 1]);
ylim(ax, [0, Matrix_Size_Acquired(3) + 1]);
xlabel(ax, 'Phase Encode (k_y)');
ylabel(ax, 'Slice Encode (k_z)');

% 3. Define Colors
unique_frames = unique(traj_table.Frame);
num_frames = length(unique_frames);
frame_colors = jet(num_frames); 

% 4. Animation Loop
for f_idx = 1:num_frames
    this_frame_num = unique_frames(f_idx);
    
    % Extract coordinates for JUST this frame
    % We find the rows for the current frame
    rows = find(traj_table.Frame == this_frame_num);
    
    % Create a NEW animated line for this specific frame
    this_color = frame_colors(f_idx, :);
    h_anim = animatedline('LineStyle', 'none', 'Marker', '.', ...
                          'MarkerSize', 12, 'Color', this_color);
    
    title(ax, sprintf('Acquiring Frame %d', this_frame_num), 'Color', this_color);
    
    % -- INTELLIGENT STRIDE SETTING --
    % Frame 0 usually has thousands of points (Prep scan). 
    % We plot it in larger chunks so it doesn't take forever.
    % Frames 1+ are usually sparse (Undersampled). We plot them slowly.
    if this_frame_num == 0
        stride = 10;   % Plot 50 points at a time for the big prep scan
        pause_time = 5e-3;
    else
        stride = 1;    % Plot 1 point at a time for the TWIST frames
        pause_time = 5e-3;  % Slow down to see the spiral pattern
    end

    % -- INNER LOOP: Animate points appearing --
    num_points = length(rows);
    for k = 1:stride:num_points
        % Determine batch range
        idx_end = min(k + stride - 1, num_points);
        batch_indices = rows(k:idx_end);
        
        % Get data
        x_data = traj_table.("Row (phase)")(batch_indices);
        y_data = traj_table.("Column (slice)")(batch_indices);
        
        % Add points
        addpoints(h_anim, x_data, y_data);
        
        % FORCE UPDATE (Key Fix)
        drawnow; 
        
        % Pause to allow eye to track motion
        if pause_time > 0
            pause(pause_time); 
        end
    end
end

title(ax, 'Acquisition Complete', 'Color', 'w');
