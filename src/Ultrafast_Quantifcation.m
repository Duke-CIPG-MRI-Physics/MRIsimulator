tic
clear; clc; close all;
load("SimulationParameters.mat")

pBs = [.1, .25, .33, .5];
pAs = .04:.1:1;

num_pBs = length(pBs);
num_pAs = length(pAs);

% 1. Pre-allocate a 2D structure array. 
% parfor requires outputs to be correctly "sliced" based on the loop index.
emptyStruct = struct('pB', [], 'pA', [], 'Measured_Contrast', [], 'Timepoints', []);
resultsStruct = repmat(emptyStruct, num_pBs, num_pAs);

% 2. Use integer indices for the parfor loop

% parpool("Threads",2);
for i = 1:num_pBs 
    % Extract the current pB value
    pB_val = pBs(i);
    
    % 3. Create a local copy of the broadcast variable (SimulationParameters)
    % This prevents workers from attempting to modify a shared struct simultaneously
    localSimParams = SimulationParameters; 
    
    for j = 1:num_pAs
        pA_val = pAs(j);
        
        localSimParams.TWIST.pB = pB_val;
        localSimParams.TWIST.pA = pA_val;
        
        Output = Analytical_TWIST_Simulator(localSimParams);
        
        % 4. Save results using the sliced indices (i, j) instead of a counter
        resultsStruct(i, j).pB = pB_val;
        resultsStruct(i, j).pA = pA_val;
        resultsStruct(i, j).Measured_Contrast = Output.measured_contrast;
        resultsStruct(i, j).Timepoints = Output.timepoints;

    end
end

% 5. Flatten the 2D struct array back into a 1D array 
% This matches the structure format of your original code
resultsStruct = resultsStruct(:);

save("test_results","resultsStruct")

toc

% Example to retrieve data later:
% To get the contrast vector for the 3rd simulation run:
% myVector = resultsStruct(3).Measured_Contrast;
%%  --- Plotting
figure
for ii = 1:height(resultsStruct)
    hold on
    plot(resultsStruct(ii).Timepoints,resultsStruct(ii).Measured_Contrast)
end
plot(1:1:1000,breastPhantomParams.lesionIntensityFunction(1:1:1000) ...
    + breastPhantomParams.breastIntensity)
hold off