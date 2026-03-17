clear; clc; %close all;
load("SimulationParameters.mat")

pBs = [0,.1,.25];
pAs = [.04:.01:.3];
n_shifts = 5;

num_pBs = length(pBs);
num_pAs = length(pAs);

% Pre-allocate with the simulated fields
emptyStruct = struct('pB', [], 'pA', [], 'delay', [], ...
    'Measured_Contrast', [], 'Timepoints', [], ...
    'Sim_Contrast', [], 'Sim_Timepoints', []);

resultsStruct = repmat(emptyStruct, num_pBs, num_pAs, n_shifts); 

SimulationParameters.TWIST.N_measurements = 15;

for i = 1:num_pBs 
    pB_val = pBs(i);
    
    for j = 1:num_pAs
        pA_val = pAs(j);
        
        localSimParams = SimulationParameters; 
        shift_amount = 0;

        for k = 1:n_shifts
            localSimParams.LesionParameters.lesionArrivalDelay_s = SimulationParameters.LesionParameters.lesionArrivalDelay_s + (shift_amount*k);
            fprintf("Simulating pA = %g, pB = %g, Delay = %g... \n", pA_val, pB_val, localSimParams.LesionParameters.lesionArrivalDelay_s)
            
            localSimParams.TWIST.pB = pB_val;
            localSimParams.TWIST.pA = pA_val;
            
            if pA_val < 0
                Output = GPU_Analytical_TWIST_Simulator(localSimParams);
            else
                Output = Analytical_TWIST_Simulator(localSimParams);
            end
            
           if k ==1 
               temporal_resolution = ...
               Output.measured.timepoints(3)-Output.measured.timepoints(2);
               shift_amount = temporal_resolution/(n_shifts-1);
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

for i = 1:num_pBs 
    pB_val = pBs(i);

    for j = 1:num_pAs
        pA_val = pAs(j);
% 1. Extract and concatenate into 1D row vectors

raw_Measured_Contrast = [resultsStruct(i, j, :).Measured_Contrast];
raw_Timepoints = [resultsStruct(i, j, :).Timepoints];

% 2. Sort the timepoints chronologically 
[combined_Timepoints, sort_idx] = sort(raw_Timepoints);

% 3. Apply the same sorting index to the contrast data
combined_Measured_Contrast = raw_Measured_Contrast(sort_idx);

            resultsStruct2(i, j).pB = pB_val;
            resultsStruct2(i, j).pA = pA_val;
            
            % Corrected assignments
            resultsStruct2(i, j).Measured_Contrast = combined_Measured_Contrast;
            resultsStruct2(i, j).Timepoints = combined_Timepoints;
                   
    end
end
%% --- Save and Flatten ---
% Save the flattened version as you originally did, but use a new variable 
% name so you don't overwrite the useful 3D struct in your workspace.
resultsStruct_flat = resultsStruct(:);
resultsStruct2_flat = resultsStruct2(:);
save("Results_test", "resultsStruct2_flat")


%% --- Plotting ---
figure;
t = tiledlayout(1, 1, 'TileSpacing', 'compact');
ax = nexttile(t);
hold(ax, 'on');

% 1. Plot the Single Simulated Ground Truth Curve FIRST (from k=1)
plot(ax, resultsStruct(1, 1, 1).Sim_Timepoints, resultsStruct(1, 1, 1).Sim_Contrast, ...
    'k--', 'LineWidth', 2, 'DisplayName', 'Analytical Ground Truth');

% 2. Plot the new combined vectors across all K
meas_name = sprintf('Combined TWIST (pA=%.2f, pB=%.2f)', pAs(1), pBs(1));
plot(ax, combined_Timepoints, combined_Measured_Contrast, ...
    '-o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', meas_name);

hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'Time Since Injection (s)', 'FontWeight', 'bold');
ylabel(ax, 'ROI Contrast Intensity', 'FontWeight', 'bold');
title(ax, 'TWIST Kinetics: Reconstructed vs. True Enhancement');
legend(ax, 'Location', 'bestoutside');