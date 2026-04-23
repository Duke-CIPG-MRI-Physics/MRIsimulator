clear; clc; close all;

% --- Inputs ---
pA = 0.04; 
pB = .1; 
Matrix_Size_Acquired = [1,256,256];
FOV_acquired = [1,256,256];
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

sampled_mask_swirl = zeros(size(sampled_mask(:,:,1)));

for i =1:size(sampled_mask,3)
    sampled_mask_swirl = sampled_mask_swirl + (double(sampled_mask(:,:,i)) * i);
end
sampled_mask_swirl(sampled_mask_swirl == max(sampled_mask_swirl)) = 0;
imshow(sampled_mask_swirl,[])

sliceViewer(double(sampled_mask)); 
title('TWIST Sampling Mask (Scroll through Frames)');

%% --- Play and Save Animation (GIF) ---
sampled_mask = sampled_mask(:,:,2:end);
nFrames = size(sampled_mask,3);
sampled_mask = uint8(sampled_mask);
sampled_mask = sampled_mask.*imread("FFT of Ultrafast.gif");

% --- Save Raw Matrix Animation (GIF) ---
gif_filename = 'TWIST_matrix_only.gif';
delay_time = 0.1; % Time in seconds between frames

if exist(gif_filename, 'file')
    delete(gif_filename);
end

% Create a simple colormap: Row 1 is Black (for 0), Row 2 is White (for 1)
cm = [0 0 0; 1 1 1]; 

for f = 1:nFrames
    % Extract the slice and convert logical to uint8 (0 and 1)
    matrix_frame = sampled_mask(:, :, f);
    
    % Write directly to the GIF File
    if f == 1
        imwrite(matrix_frame, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(matrix_frame, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end
disp(['Raw matrix GIF saved successfully as: ', gif_filename]);