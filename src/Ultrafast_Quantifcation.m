tic
clear; clc; close all;
load("SimulationParameters.mat")

pBs = [0];
pAs = .04:.01:1;

num_pBs = length(pBs);
num_pAs = length(pAs);

% Pre-allocate with the simulated fields
emptyStruct = struct('pB', [], 'pA', [], ...
    'Measured_Contrast', [], 'Timepoints', [], ...
    'Sim_Contrast', [], 'Sim_Timepoints', []);
resultsStruct = repmat(emptyStruct, num_pBs, num_pAs);

for i = 1:num_pBs 
    pB_val = pBs(i);
    localSimParams = SimulationParameters; 
    
    for j = 1:num_pAs
        pA_val = pAs(j);
        fprintf("Simulating pA = %g, pB = %g...\n",pA_val,pB_val)
        localSimParams.TWIST.pB = pB_val;
        localSimParams.TWIST.pA = pA_val;

        if pA < .6
        Output = GPU_Analytical_TWIST_Simulator(localSimParams);
        else
        Output = Analytical_TWIST_Simulator(localSimParams);
        end

        
        resultsStruct(i, j).pB = pB_val;
        resultsStruct(i, j).pA = pA_val;
        
        % Corrected assignments
        resultsStruct(i, j).Measured_Contrast = Output.measured.contrast;
        resultsStruct(i, j).Timepoints = Output.measured.timepoints;
        
        % Capture the shifting simulated ground truth
        resultsStruct(i, j).Sim_Contrast = Output.simulated.contrast;
        resultsStruct(i, j).Sim_Timepoints = Output.simulated.timepoints;
        fprintf("Done\n")
    end
end
resultsStruct = resultsStruct(:);

save("Results_pB_0","resultsStruct")

toc

% Example to retrieve data later:
% To get the contrast vector for the 3rd simulation run:
% myVector = resultsStruct(3).Measured_Contrast;

%%  --- Plotting
figure;
t = tiledlayout(1, 1, 'TileSpacing', 'compact');
ax = nexttile(t);
hold(ax, 'on');

% Plot the Single Simulated Ground Truth Curve FIRST
plot(ax, resultsStruct(1).Sim_Timepoints, resultsStruct(1).Sim_Contrast, ...
    'k--', 'LineWidth', 2, 'DisplayName', 'Analytical Ground Truth');

% Define a color order for the measured curves
colors = lines(length(resultsStruct));

% Loop through and overlay the measured TWIST data
for ii = 1:length(resultsStruct)
    meas_name = sprintf('TWIST (pA=%.2f, pB=%.2f)', resultsStruct(ii).pA, resultsStruct(ii).pB);
    
    plot(ax, resultsStruct(ii).Timepoints, resultsStruct(ii).Measured_Contrast, ...
        '-o', 'Color', colors(ii,:), 'LineWidth', 1.5, 'DisplayName', meas_name);
end

hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'Time Since Injection (s)', 'FontWeight', 'bold');
ylabel(ax, 'ROI Contrast Intensity', 'FontWeight', 'bold');
title(ax, 'TWIST Kinetics: Reconstructed vs. True Enhancement');
legend(ax, 'Location', 'bestoutside');