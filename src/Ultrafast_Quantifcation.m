clear; clc; close all;
load("SimulationParameters.mat")

pBs = [.1,.25];
pAs = [.04:.01:.15];
n_shifts = 5;

num_pBs = length(pBs);
num_pAs = length(pAs);

% Pre-allocate with the simulated fields
emptyStruct = struct('pB', [], 'pA', [], 'delay', [], ...
    'Measured_Contrast_L', [], 'Measured_Contrast_M', [], ...
    'Measured_Contrast_S', [], 'Measured_Contrast_XS', [], ...
    'Timepoints', [], 'Sim_Contrast', [], 'Sim_Timepoints', []);

resultsStruct = repmat(emptyStruct, num_pBs, num_pAs, n_shifts); 

% SimulationParameters.TWIST.N_measurements = 15;

for i = 1:num_pBs 
    pB_val = pBs(i);
    
    for j = 1:num_pAs
        pA_val = pAs(j);
        
        localSimParams = SimulationParameters; 
        shift_amount = 0;

        for k = 1:n_shifts
            localSimParams.LesionParameters.lesionArrivalDelay_s = ...
            SimulationParameters.LesionParameters.lesionArrivalDelay_s + (shift_amount * (k-1));
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
            resultsStruct(i, j, k).Measured_Contrast_L = Output.measured.contrast_L;
            resultsStruct(i, j, k).Measured_Contrast_M = Output.measured.contrast_M;
            resultsStruct(i, j, k).Measured_Contrast_S = Output.measured.contrast_S;
            resultsStruct(i, j, k).Measured_Contrast_XS = Output.measured.contrast_XS;
            resultsStruct(i, j, k).Timepoints = Output.measured.timepoints;
            resultsStruct(i, j, k).nominal_temporal_resolution = Output.measured.nominal_temporal_resolution;

            % Capture the shifting simulated ground truth
            resultsStruct(i, j, k).Sim_Contrast = Output.simulated.contrast;
            resultsStruct(i, j, k).Sim_Timepoints = Output.simulated.timepoints;
            fprintf("Done\n")
        end
    end
end

%% Extract Data Across K for a specific I and J ---
resultsStruct2 = struct('pB', [], 'pA', [], 'delay', [], ...
    'Measured_Contrast_L', [], 'Measured_Contrast_M', [], ...
    'Measured_Contrast_S', [], 'Measured_Contrast_XS', [], ...
    'Timepoints', [], 'Sim_Contrast', [], 'Sim_Timepoints', []);

for i = 1:num_pBs 
    pB_val = pBs(i);

    for j = 1:num_pAs
        pA_val = pAs(j);
% 1. Extract and concatenate into 1D row vectors
raw_Measured_Contrast_L = [resultsStruct(i, j, :).Measured_Contrast_L];
raw_Measured_Contrast_M = [resultsStruct(i, j, :).Measured_Contrast_M];
raw_Measured_Contrast_S = [resultsStruct(i, j, :).Measured_Contrast_S];
raw_Measured_Contrast_XS = [resultsStruct(i, j, :).Measured_Contrast_XS];

raw_Timepoints = [resultsStruct(i, j, :).Timepoints];

% 2. Sort the timepoints chronologically 
[combined_Timepoints, sort_idx] = sort(raw_Timepoints);

% 3. Apply the same sorting index to the contrast data
combined_Measured_Contrast_L = raw_Measured_Contrast_L(sort_idx);
combined_Measured_Contrast_M = raw_Measured_Contrast_M(sort_idx);
combined_Measured_Contrast_S = raw_Measured_Contrast_S(sort_idx);
combined_Measured_Contrast_XS = raw_Measured_Contrast_XS(sort_idx);



            resultsStruct2(i, j).pB = pB_val;
            resultsStruct2(i, j).pA = pA_val;
            
            % Corrected assignments
            resultsStruct2(i, j).Measured_Contrast_L = combined_Measured_Contrast_L;
            resultsStruct2(i, j).Measured_Contrast_M = combined_Measured_Contrast_M;
            resultsStruct2(i, j).Measured_Contrast_S = combined_Measured_Contrast_S;
            resultsStruct2(i, j).Measured_Contrast_XS = combined_Measured_Contrast_XS;
            resultsStruct2(i, j).Timepoints = combined_Timepoints;
            resultsStruct(i,j).nominal_temporal_resolution = resultsStruct(i,j).nominal_temporal_resolution;
                   
    end
end

%% --- Save and Flatten ---
% Save the flattened version as you originally did, but use a new variable 
% name so you don't overwrite the useful 3D struct in your workspace.
resultsStruct_flat = resultsStruct(:);
resultsStruct2_flat = resultsStruct2(:);
Simulated.contrast = Output.simulated.contrast;
Simulated.timepoints = Output.simulated.timepoints;
save("Results_forward", "resultsStruct2_flat","Simulated")

