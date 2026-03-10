tic
clear; clc; %close all;
load("SimulationParameters.mat")

pBs = [0,.1];
pAs = [.04];

num_pBs = length(pBs);
num_pAs = length(pAs);

% Pre-allocate with the simulated fields
emptyStruct = struct('pB', [], 'pA', [], ...
    'Measured_Contrast', [], 'Timepoints', [], ...
    'Sim_Contrast', [], 'Sim_Timepoints', []);

% Note: We pre-allocate to 3D since your k loop goes 1:2
resultsStruct = repmat(emptyStruct, num_pBs, num_pAs, 2); 

SimulationParameters.TWIST.N_measurements = 10;

for i = 1:num_pBs 
    pB_val = pBs(i);
    
    for j = 1:num_pAs
        pA_val = pAs(j);
        
        localSimParams = SimulationParameters; 
       
        for k = 1:2
            localSimParams.LesionParameters.lesionArrivalDelay_s = localSimParams.LesionParameters.lesionArrivalDelay_s + 5;
            fprintf("Simulating pA = %g, pB = %g, Delay = %g... \n", pA_val, pB_val, localSimParams.LesionParameters.lesionArrivalDelay_s)
            
            localSimParams.TWIST.pB = pB_val;
            localSimParams.TWIST.pA = pA_val;
            
            if pA_val < 0
                Output = GPU_Analytical_TWIST_Simulator(localSimParams);
            else
                Output = Analytical_TWIST_Simulator(localSimParams);
            end
            
            resultsStruct(i, j, k).pB = pB_val;
            resultsStruct(i, j, k).pA = pA_val;
            resultsStruct(i, j, k).delay = localSimParams.LesionParameters.lesionArrivalDelay_s;
            
            % Corrected assignments
            resultsStruct(i, j, k).Measured_Contrast = Output.measured.contrast;
            resultsStruct(i, j, k).Timepoints = Output.measured.timepoints;
                   
            % Capture the shifting simulated ground truth
            resultsStruct(i, j, k).Sim_Contrast = Output.simulated.contrast;
            resultsStruct(i, j, k).Sim_Timepoints = Output.simulated.timepoints;
            fprintf("Done\n")
        end
    end
end

%% Extract Data Across K for a specific I and J ---
% Specify the target combination you want to plot
target_i = 1; %TODO: make this loop
target_j = 1; 

% 1. Extract and concatenate into 1D row vectors
raw_Measured_Contrast = [resultsStruct(target_i, target_j, :).Measured_Contrast];
raw_Timepoints = [resultsStruct(target_i, target_j, :).Timepoints];

% 2. Sort the timepoints chronologically so the plot line doesn't zig-zag
[combined_Timepoints, sort_idx] = sort(raw_Timepoints);

% 3. Apply the same sorting index to the contrast data
combined_Measured_Contrast = raw_Measured_Contrast(sort_idx);

%% --- Save and Flatten ---
% Save the flattened version as you originally did, but use a new variable 
% name so you don't overwrite the useful 3D struct in your workspace.
resultsStruct_flat = resultsStruct(:);
save("Results_test", "resultsStruct_flat")
toc

%% --- Plotting ---
figure;
t = tiledlayout(1, 1, 'TileSpacing', 'compact');
ax = nexttile(t);
hold(ax, 'on');

% 1. Plot the Single Simulated Ground Truth Curve FIRST (from k=1)
plot(ax, resultsStruct(target_i, target_j, 1).Sim_Timepoints, resultsStruct(target_i, target_j, 1).Sim_Contrast, ...
    'k--', 'LineWidth', 2, 'DisplayName', 'Analytical Ground Truth');

% 2. Plot the new combined vectors across all K
meas_name = sprintf('Combined TWIST (pA=%.2f, pB=%.2f)', pAs(target_j), pBs(target_i));
plot(ax, combined_Timepoints, combined_Measured_Contrast, ...
    '-o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', meas_name);

hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'Time Since Injection (s)', 'FontWeight', 'bold');
ylabel(ax, 'ROI Contrast Intensity', 'FontWeight', 'bold');
title(ax, 'TWIST Kinetics: Reconstructed vs. True Enhancement');
legend(ax, 'Location', 'bestoutside');