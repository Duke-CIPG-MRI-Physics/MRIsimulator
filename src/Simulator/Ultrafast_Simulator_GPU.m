clear;clc;close all

%% --- Ground-Truth Segmentation and MRI 
%convention for dimensions is: [frequency,phase,slice,timepoint,coil]
fprintf('Getting Segmentation and Simulating GRE...')
Type = 3;

flip_angle_degrees = 20;
TE_s = 2.46e-3;
TR_s = 6e-3;

xcat_MRI_real = GRE_Simulator_XCAT(flip_angle_degrees,TE_s,TR_s);
%% --- Downsampling
fprintf('\nDownsampling...')
% Define desired final matrix size
Desired_Matrix_Size = [256, 256, 256]; % [frequency, phase, slice]

% 1. Perform 3D FFT on ALL timepoints at once.
% This is done by applying 1D FFTs and shifts sequentially along the first 3 dimensions.
xcat_MRI_real = gpuArray(xcat_MRI_real);
k_temp = fft(fft(fft(xcat_MRI_real, [], 1), [], 2), [], 3);
xcat_MRI_kspace = fftshift(fftshift(fftshift(k_temp, 1), 2), 3);
[xcat_MRI_kspace,xcat_MRI_real,k_temp] = gather(xcat_MRI_kspace,xcat_MRI_real,k_temp);

% 2. Calculate the cropping ranges once.
original_kspace_size = size(xcat_MRI_kspace);

% A robust way to calculate centered crop ranges
start_indices = floor((original_kspace_size(1:3) - Desired_Matrix_Size) / 2) + 1;
end_indices = start_indices + Desired_Matrix_Size - 1;

freq_range = start_indices(1):end_indices(1);
phase_range = start_indices(2):end_indices(2);
slice_range = start_indices(3):end_indices(3);

% 3. Downsample k-space by cropping all timepoints in a single operation.
% The ':' selects all elements along the 4th dimension (time).
xcat_MRI_kspace_downsampled = xcat_MRI_kspace(freq_range, phase_range, slice_range, :);

% 4. Perform inverse 3D FFT on ALL cropped timepoints at once.
xcat_MRI_kspace_downsampled = gpuArray(xcat_MRI_kspace_downsampled);
k_shifted_back = ifftshift(ifftshift(ifftshift(xcat_MRI_kspace_downsampled, 1), 2), 3);
xcat_MRI_real_downsampled = ifft(ifft(ifft(k_shifted_back, [], 1), [], 2), [], 3);
[xcat_MRI_real_downsampled,xcat_MRI_kspace_downsampled] = gather(xcat_MRI_real_downsampled,xcat_MRI_kspace_downsampled);

% 5. Final assignment.
xcat_MRI_real = abs(xcat_MRI_real_downsampled);
xcat_MRI_kspace = xcat_MRI_kspace_downsampled;

clear('k_temp')
%% --- Simulating Coils and adding noise
%Eventually we'll want to just load shading functions, no need to calculate
%every time
fprintf('\nSimulating Coils...')
nCoilsEven = 4;
sigma = [0.8 0.55];
coilSplitDirection = 'z';
coilDirection = 'y'; % This might have to change
opts = struct;
opts.planeOffset = 1.02;
opts.phaseAmp = 0.6;
opts.noiseRMS = 1e-2;     % adjust for desired SNR
opts.dtype = "single";      % or "double"    (must be *string*)

if Type == 3
   xcat_coils_real = zeros([size(xcat_MRI_real),1,nCoilsEven]);
   xcat_coils_kspace = xcat_coils_real;
elseif Type == 4
    xcat_coils_real = zeros([size(xcat_MRI_real),nCoilsEven]);
    xcat_coils_kspace = xcat_coils_real;
end


for ii_timepoint = 1:size(xcat_MRI_real,4)
    [xcat_coils_real(:,:,:,ii_timepoint,:),xcat_coils_kspace(:,:,:,ii_timepoint,:),~] = ...
    simulateCoils(xcat_MRI_real(:,:,:,ii_timepoint),nCoilsEven,sigma,coilSplitDirection,coilDirection);
end
%% --- TWIST Point Identification and Ordering
fprintf('\nIdentifying TWIST Points...\n')
pause(.2)
pA = .3;
N = 10;
Time_Measured = 1;
TR = 6e-3;
R = [2,2];
PF_Factor = [1,1];
[Sampling_Table,PF_ask,GRAPPA_ask] = Ultrafast_Sampling(Desired_Matrix_Size,pA,N,Time_Measured,TR,R,PF_Factor);
%% --- Sampling Input Data in accordance with Sampling_Table
fprintf('\nSampling TWIST Points...')

%Here we are splitting up the sampling table based on each Bj and using the
%index values to read into the complete kspace segmentation

TWIST_kspace = zeros([Desired_Matrix_Size,max(Sampling_Table.Bj)+1,nCoilsEven]);


% Get the dimensions of your matrices once, outside the loops
[n_freq, n_phase, n_slice, n_timepoints, n_coils] = size(TWIST_kspace);
source_dims = size(xcat_coils_kspace);

% Loop over coils and timepoints
for jj_coil = 1:n_coils
    for ii_timepoint = 1:n_timepoints

        % 1. Find all table entries for the current timepoint
        is_timepoint = (Sampling_Table.Bj == (ii_timepoint - 1));

        % 2. Extract the corresponding row and column indices
        % Note: Using .("Name") syntax for column names with spaces/parentheses
        rows = Sampling_Table.("Row (phase)")(is_timepoint);
        cols = Sampling_Table.("Column (slice)")(is_timepoint);

        % If there are no samples for this timepoint, skip to the next iteration
        if isempty(rows)
            continue;
        end

        % 3. Generate expanded index vectors for all frequency points
        % This creates a full list of every source/destination element to be copied
        num_samples = numel(rows);
        freq_idx_col = (1:n_freq)'; % Creates a column vector [1; 2; ...; n_freq]

        % Create full index vectors for each dimension
        freq_full = repmat(freq_idx_col, num_samples, 1);
        rows_full = repelem(rows(:), n_freq, 1);
        cols_full = repelem(cols(:), n_freq, 1);

        % 4. Convert the subscript indices to linear indices for both matrices
        % This is the key step that allows for a single, vectorized assignment

        % Indices for the source matrix (xcat_coils_kspace)
        source_lin_idx = sub2ind(source_dims, ...
                                 freq_full, rows_full, cols_full, ...
                                 1, jj_coil);

        % Indices for the destination matrix (TWIST_kspace)
        dest_lin_idx = sub2ind(size(TWIST_kspace), ...
                                freq_full, rows_full, cols_full, ...
                                ii_timepoint, jj_coil);

        % 5. Perform the assignment in a single, vectorized operation
        TWIST_kspace(dest_lin_idx) = xcat_coils_kspace(source_lin_idx);

    end
end
clear('cols_full','rows_full','freq_full','source_lin_idx','dest_lin_idx')
%% ---  Updating K-Space from each measurement to undo TWIST
fprintf('\nUndoing TWIST...')

% Initialize the output array
unTWISTed_kspace = zeros(size(TWIST_kspace));

% The first timepoint is the baseline for all coils
unTWISTed_kspace(:,:,:,1,:) = TWIST_kspace(:,:,:,1,:);

% Loop only through timepoints (vectorized over coils)
for ii_timepoint = 2:size(TWIST_kspace,4)
    % 1. Create a temporary variable for the current timepoint's slice,
    % starting with data from the previous timepoint.
    dest_slice = unTWISTed_kspace(:,:,:,ii_timepoint-1,:);

    % 2. Get the new sparse measurements for the current timepoint
    source_slice = TWIST_kspace(:,:,:,ii_timepoint,:);

    % 3. Create a logical mask of where the new measurements exist
    update_mask = (source_slice ~= 0);

    % 4. Use the mask to update the destination slice with the new values.
    % This logical indexing is very fast.
    dest_slice(update_mask) = source_slice(update_mask);

    % 5. Assign the updated slice back into the main array.
    unTWISTed_kspace(:,:,:,ii_timepoint,:) = dest_slice;
end

clear('dest_slice','source_slice','update_mask')

% Perform the inverse 3D FFT on all timepoints and coils at once.
% This part of the optimization remains correct and highly efficient.
unTWISTed_kspace = gpuArray(unTWISTed_kspace);
k_shifted_back = ifftshift(ifftshift(ifftshift(unTWISTed_kspace, 1), 2), 3);
unTWISTed_real = ifft(ifft(ifft(k_shifted_back, [], 1), [], 2), [], 3);
[unTWISTed_kspace,unTWISTed_real] = gather(unTWISTed_kspace,unTWISTed_real);
%% --- Undoing GRAPPA
if GRAPPA_ask == 'y'
fprintf('\nUndoing GRAPPA...')
R = [1,R];

%first we need to establish a calibration region
n_calib_phases = 24;
n_calib_slices = 24;

center_phases = round((size(unTWISTed_kspace,2)-n_calib_phases)/2:((size(unTWISTed_kspace,2)-n_calib_phases)/2)+n_calib_phases-1);
center_slices = round((size(unTWISTed_kspace,3)-n_calib_slices)/2:((size(unTWISTed_kspace,3)-n_calib_slices)/2)+n_calib_slices-1);
 
% Calibration region is defined from original coil images as they are FULLY
% sampled
calibration_region = squeeze(xcat_coils_kspace(:,center_phases,center_slices,1,:));
calibration_region = permute(calibration_region,[4,1,2,3]);

Output_kspace = zeros(size(unTWISTed_kspace,5),size(unTWISTed_kspace,1),size(unTWISTed_kspace,2),size(unTWISTed_kspace,3),size(unTWISTed_kspace,4));

for ii_timepoint = 1:size(unTWISTed_kspace,4)
    %need to permute because grappa expect coils as first dim
    input = squeeze(unTWISTed_kspace(:,:,:,ii_timepoint,:));
    input = permute(input,[4,1,2,3]);
         
    Output_kspace(:,:,:,:,ii_timepoint)   = grappa(input,calibration_region, R , [3,3,3]);
end

Output_Coils_K = permute(Output_kspace,[2,3,4,5,1]);

else
Output_Coils_K = unTWISTed_kspace;
end
clear('Output_kspace','input','calibration_region')
%% --- IFFT and Sum of Squares
fprintf('\nPerforming IFFT and Sum of Squares...')

% Perform the inverse 3D FFT on all timepoints and coils simultaneously
% by applying 1D functions sequentially along the first 3 dimensions.
% This replaces the first set of nested loops.
Output_Coils_K = gpuArray(Output_Coils_K);
k_shifted_back = ifftshift(ifftshift(ifftshift(Output_Coils_K, 1), 2), 3);
Output_Coils_Real = ifft(ifft(ifft(k_shifted_back, [], 1), [], 2), [], 3);
[Output_Coils_K,Output_Coils_Real] = gather(Output_Coils_K,Output_Coils_Real);

% Calculate the root sum of squares directly on the 5D array.
% The sum() function operates along the specified dimension (5, for coils),
% replacing the second loop. The 'abs()' ensures the magnitude of the
% complex signal is used.
SOS = sqrt(sum(abs(Output_Coils_Real).^2, 5));
clear('k_shifted_back')
%% ---  Displaying
fprintf('\nDisplaying...')
figure
imshow3D(log(abs(squeeze(TWIST_kspace(50,:,:,:,1)))))
title('K-space being acquired each measurement, Coil 1')

figure
imshow3D(log(abs(squeeze(unTWISTed_kspace(50,:,:,:,1)))))
title('K-space being updated with each measurement, Coil 1')

figure
imshow3D(abs(squeeze(unTWISTed_real(50,:,:,:,1))))

figure
imshow3D(log(abs(squeeze(Output_Coils_K(50,:,:,:,1)))))

imslice(abs(Output_Coils_Real))
title('Coil Images')

imslice(abs(SOS))
title('Sum of Squares')






