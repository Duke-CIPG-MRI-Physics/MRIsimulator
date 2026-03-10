clear; clc; close all;

% --- Inputs ---
pA = 0.05; 
pB = 0.1; 
Matrix_Size_Acquired = [1,100,100];
FOV_acquired = [1,200,100];
R = [1,1];
PF_Factor = [1,1]; 

% --- Run TWIST ---
[TWIST_sampling_order] = TWIST(pA, pB, Matrix_Size_Acquired, FOV_acquired, R, PF_Factor);

% --- Correct Visualization Setup ---
nRows = Matrix_Size_Acquired(2);
nCols = Matrix_Size_Acquired(3);
nFrames = max(TWIST_sampling_order.Frame) + 1; % +1 because frames start at 0

% Initialize 3D mask: [Phase, Slice, Time]
sampled_mask = false(nRows, nCols, nFrames);

% Loop through each row of the table to place points in the correct "Time" slice
for i = 1:height(TWIST_sampling_order)
    r = TWIST_sampling_order.("Row (phase)")(i);
    c = TWIST_sampling_order.("Column (slice)")(i);
    f = TWIST_sampling_order.Frame(i) + 1; % Shift 0-index to 1-index for MATLAB
    
    sampled_mask(r, c, f) = true;
end


sliceViewer(double(sampled_mask)); 
title('TWIST Sampling Mask (Scroll through Frames)');