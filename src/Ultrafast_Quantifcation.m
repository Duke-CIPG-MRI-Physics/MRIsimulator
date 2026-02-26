%This script is designed for performance testing on simulated
%Ultrafast breast MRI

clear; clc; close all;
load("SimulationParameters.mat")

pBs = [.1, .25, .33, .5];
pAs = .04:.2:1;

% Pre-allocate a 2D cell array (rows = pBs, cols = pAs)
resultsCell = cell(length(pBs), length(pAs));

% Loop using indices instead of values
for idx_B = 1:length(pBs)
    for idx_A = 1:length(pAs)
        
        % Extract the actual values for the simulation
        current_pB = pBs(idx_B);
        current_pA = pAs(idx_A);
        
        SimulationParameters.TWIST.pB = current_pB;
        SimulationParameters.TWIST.pA = current_pA;
        Output = Analytical_TWIST_Simulator(SimulationParameters);
        
        % Store the resulting vector in the corresponding grid cell
        resultsCell{idx_B, idx_A} = Output.measured_contrast;
    end
end

% Example to retrieve data later:
% To get the contrast vector where pB is .25 (index 2) and pA is .44 (index 3):
% myVector = resultsCell{2, 3};


%% - Deviation

percent_deviation_from_GT = 100*abs((contrast_values_measured-gt_at_TWIST_pts)./gt_at_TWIST_pts);

max_percent_deviation_from_GT = max(percent_deviation_from_GT);
fprintf("Max deviation from ground-truth: %g%%\n",max_percent_deviation_from_GT)
average_percent_deviation_from_GT = mean(percent_deviation_from_GT);
fprintf("Avg. deviation from ground-truth: %g%%\n",average_percent_deviation_from_GT)
