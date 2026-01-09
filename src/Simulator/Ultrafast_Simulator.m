clear;clc;close all
% Roberto Carrascosa
% Duke University, Medical Phyiscs
% 2025


%Below is a high-level overview of each section of this script

% 1. Ground-Truth Segmentation and Bloch Simulation:
%       This yields a 3D or 4D MRI volume based on Bloch Simulations. It is
%       a 'perfect' MRI which we will reference later when we need to
%       sample points of k-space
%--------------------------------------------------------
% 2. Downsampling:
%       The variables get huge for large input matrices so we provide the
%       option to downsample to any arbitrary matrix size for ease of use
%       and efficiency when debugging. Works for both 3D and 4D, but
%       obviously doesn't downsample across time dimension
%--------------------------------------------------------
% 3. Simulating Coils and adding Noise
%       An actual scanner doesn't capture the entire volume at once, there
%       are several coils whose images are combined to create a single
%       volume. This function simulates those individual coil elements as
%       they are needed for later steps.
%
%       This also adds noise to create a more realistic volume
%--------------------------------------------------------
% 4. TWIST Point Identification and Ordering, GRAPPA and PF Undersampling
%       This part doesn't reference the image volume at all. It simply uses
%       the given matrix size to generate an ordered list of
%       points/coordinates in the matrix for sampling according to TWIST
%
%       CRUCIALLY, this is also the section where we remove lines
%       according to GRAPPA and Partial Fourier
%--------------------------------------------------------
% 5. Sampling Input Data in accordance with Sampling_Table
%       Here we actually sample the image volume according to the ordered
%       list we created in the previous section. It separates them
%       according to the TWIST time points
%
%       If GRAPPA or PF was called for in previos section, this output
%       volume will be undersampled accordingly.
%--------------------------------------------------------
%--------------------------------------------------------

% CRITICAL SHIFT: Everything up to this point has been about generating a
% 'realistic' MRI that has been undersampled and divided into time points
% according to TWIST, GRAPPA, and PF. From here forward, we will perform
% reconstruction to undo those undersampling techniques. 

%--------------------------------------------------------
%--------------------------------------------------------
% 6. Updating K-Space from each measurement to undo TWIST
%       This section does the sliding window reconstruction to 'undo' TWIST
%       and fill-in k-space
%--------------------------------------------------------
% 7. Undoing GRAPPA
%       Pretty self-explanatory, this is the section which needs the coil
%       images.
%--------------------------------------------------------
% 8. Undoing Partial  Fourier
%
%--------------------------------------------------------
% 9.IFFT and Sum of Squares
%       This section combines the un-GRAPPAd coil images and goes back to
%       image space



%% --- 1. Ground-Truth Segmentation and Bloch Simulation
%convention for dimensions is: [frequency,phase,slice,timepoint,coil]
fprintf('Getting Segmentation and Simulating GRE...')
Type = '3D';

if all(Type == '3D')
%This simulates a single 3D volume

    flip_angle_degrees = 20;
    TE_s = 2.46e-3;
    TR_s = 6e-3;

    xcat_MRI_IMspace = GRE_Simulator_XCAT(flip_angle_degrees,TE_s,TR_s);

elseif all(Type == '4D')
%This loads a premade 4D volume

    load("xcat_MRI_4D.mat")
    xcat_MRI_IMspace = xcat_MRI_4D;

else
error("Type must be '3D' or '4D'")
end
   


%% --- 2. Downsampling
fprintf('\nDownsampling...')
% Define desired final matrix size
Desired_Matrix_Size = [128, 128, 128]; % [frequency, phase, slice]

% 1. Perform 3D FFT on ALL timepoints at once.
% This is done by applying 1D FFTs and shifts sequentially along the first 3 dimensions.
k_temp = fft(fft(fft(xcat_MRI_IMspace, [], 1), [], 2), [], 3);
xcat_MRI_Kspace = fftshift(fftshift(fftshift(k_temp, 1), 2), 3);

% 2. Calculate the cropping ranges once.
original_Kspace_size = size(xcat_MRI_Kspace);

% A robust way to calculate centered crop ranges
start_indices = floor((original_Kspace_size(1:3) - Desired_Matrix_Size) / 2) + 1;
end_indices = start_indices + Desired_Matrix_Size - 1;

freq_range = start_indices(1):end_indices(1);
phase_range = start_indices(2):end_indices(2);
slice_range = start_indices(3):end_indices(3);

% 3. Downsample k-space by cropping all timepoints in a single operation.
% The ':' selects all elements along the 4th dimension (time).
xcat_MRI_Kspace_downsampled = xcat_MRI_Kspace(freq_range, phase_range, slice_range, :);

% 4. Perform inverse 3D FFT on ALL cropped timepoints at once.
k_shifted_back = ifftshift(ifftshift(ifftshift(xcat_MRI_Kspace_downsampled, 1), 2), 3);
xcat_MRI_IMspace_downsampled = ifft(ifft(ifft(k_shifted_back, [], 1), [], 2), [], 3);

% 5. Final assignment.
xcat_MRI_IMspace = abs(xcat_MRI_IMspace_downsampled);
xcat_MRI_Kspace = xcat_MRI_Kspace_downsampled;

clear('k_temp','xcat_MRI_IMspace_downsampled','xcat_MRI_Kspace_downsampled','k_shifted_back')
%% --- 3. Simulating Coils and adding Noise
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

if all(Type == '3D')
   xcat_coils_IMspace = zeros([size(xcat_MRI_IMspace),1,nCoilsEven]);
   xcat_coils_Kspace = xcat_coils_IMspace;
elseif all(Type == '4D')
    xcat_coils_IMspace = zeros([size(xcat_MRI_IMspace),nCoilsEven]);
    xcat_coils_Kspace = xcat_coils_IMspace;
end

%This loop calls the function to simulate coils
for ii_timepoint = 1:size(xcat_MRI_IMspace,4)
    [xcat_coils_IMspace(:,:,:,ii_timepoint,:),xcat_coils_Kspace(:,:,:,ii_timepoint,:),~] = ...
    simulateCoils(xcat_MRI_IMspace(:,:,:,ii_timepoint),nCoilsEven,sigma,coilSplitDirection,coilDirection);
end
%% --- 4. TWIST Point Identification and Ordering, GRAPPA and PF Undersampling
fprintf('\nIdentifying TWIST Points...\n')
pause(.2)

pA = .1;
N = 10;
Time_Measured = 1;
TR = 6e-3;
R = [2,2];
PF_Factor = [.8,.8];

[Sampling_Table,PF_ask,GRAPPA_ask] = Ultrafast_Sampling(Desired_Matrix_Size,pA,N,Time_Measured,TR,R,PF_Factor);
%% --- 5. Sampling Input Data in accordance with Sampling_Table
fprintf('\nSampling TWIST Points...')

%Here we are splitting up the sampling table based on each Bj and using the
%index values to read into the complete kspace segmentation

TWIST_Kspace = zeros([Desired_Matrix_Size,max(Sampling_Table.Bj)+1,nCoilsEven]);


% Get the dimensions of your matrices once, outside the loops
[n_freq, n_phase, n_slice, n_timepoints, n_coils] = size(TWIST_Kspace);
source_dims = size(xcat_coils_Kspace);

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

        % Indices for the source matrix (xcat_coils_Kspace)
        source_lin_idx = sub2ind(source_dims, ...
                                 freq_full, rows_full, cols_full, ...
                                 1, jj_coil);

        % Indices for the destination matrix (TWIST_Kspace)
        dest_lin_idx = sub2ind(size(TWIST_Kspace), ...
                                freq_full, rows_full, cols_full, ...
                                ii_timepoint, jj_coil);

        % 5. Perform the assignment in a single, vectorized operation
        TWIST_Kspace(dest_lin_idx) = xcat_coils_Kspace(source_lin_idx);

    end
end
clear('cols_full','rows_full','freq_full','source_lin_idx','dest_lin_idx')
%% ---  6. Updating K-Space from each measurement to undo TWIST
fprintf('\nUndoing TWIST...')

% Initialize the output array
unTWISTed_Kspace = zeros(size(TWIST_Kspace));

% The first timepoint is the baseline for all coils
unTWISTed_Kspace(:,:,:,1,:) = TWIST_Kspace(:,:,:,1,:);

% Loop only through timepoints (vectorized over coils)
for ii_timepoint = 2:size(TWIST_Kspace,4)
    % 1. Create a temporary variable for the current timepoint's slice,
    % starting with data from the previous timepoint.
    dest_slice = unTWISTed_Kspace(:,:,:,ii_timepoint-1,:);

    % 2. Get the new sparse measurements for the current timepoint
    source_slice = TWIST_Kspace(:,:,:,ii_timepoint,:);

    % 3. Create a logical mask of where the new measurements exist
    update_mask = (source_slice ~= 0);

    % 4. Use the mask to update the destination slice with the new values.
    % This logical indexing is very fast.
    dest_slice(update_mask) = source_slice(update_mask);

    % 5. Assign the updated slice back into the main array.
    unTWISTed_Kspace(:,:,:,ii_timepoint,:) = dest_slice;
end

clear('dest_slice','source_slice','update_mask')

% Perform the inverse 3D FFT on all timepoints and coils at once.
% This part of the optimization remains correct and highly efficient.
k_shifted_back = ifftshift(ifftshift(ifftshift(unTWISTed_Kspace, 1), 2), 3);
unTWISTed_IMspace = ifft(ifft(ifft(k_shifted_back, [], 1), [], 2), [], 3);
%% --- 7. Undoing GRAPPA
if GRAPPA_ask == 'y'
fprintf('\nUndoing GRAPPA...')
R = [1,R];

%first we need to establish a calibration region
n_calib_phases = 24;
n_calib_slices = 24;

center_phases = round((size(unTWISTed_Kspace,2)-n_calib_phases)/2:((size(unTWISTed_Kspace,2)-n_calib_phases)/2)+n_calib_phases-1);
center_slices = round((size(unTWISTed_Kspace,3)-n_calib_slices)/2:((size(unTWISTed_Kspace,3)-n_calib_slices)/2)+n_calib_slices-1);
 
% Calibration region is defined from original coil images as they are FULLY
% sampled
calibration_region = squeeze(xcat_coils_Kspace(:,center_phases,center_slices,1,:));
calibration_region = permute(calibration_region,[4,1,2,3]);

Output_Kspace = zeros(size(unTWISTed_Kspace,5),size(unTWISTed_Kspace,1),size(unTWISTed_Kspace,2),size(unTWISTed_Kspace,3),size(unTWISTed_Kspace,4));

for ii_timepoint = 1:size(unTWISTed_Kspace,4)
    %need to permute because grappa expect coils as first dim
    input = squeeze(unTWISTed_Kspace(:,:,:,ii_timepoint,:));
    input = permute(input,[4,1,2,3]);
         
    Output_Kspace(:,:,:,:,ii_timepoint)   = grappa(input,calibration_region, R , [3,3,3]);
end

Output_Coils_K = permute(Output_Kspace,[2,3,4,5,1]);

else
Output_Coils_K = unTWISTed_Kspace;
end
clear('Output_Kspace','input','calibration_region')
%% --- 9. IFFT and Sum of Squares
fprintf('\nPerforming IFFT and Sum of Squares...')

% Perform the inverse 3D FFT on all timepoints and coils simultaneously
% by applying 1D functions sequentially along the first 3 dimensions.
% This replaces the first set of nested loops.
k_shifted_back = ifftshift(ifftshift(ifftshift(Output_Coils_K, 1), 2), 3);
Output_Coils_IMspace = ifft(ifft(ifft(k_shifted_back, [], 1), [], 2), [], 3);

% Calculate the root sum of squares directly on the 5D array.
% The sum() function operates along the specified dimension (5, for coils),
% replacing the second loop. The 'abs()' ensures the magnitude of the
% complex signal is used.
SOS = sqrt(sum(abs(Output_Coils_IMspace).^2, 5));
clear('k_shifted_back')
%% ---  Displaying
fprintf('\nDisplaying...')
figure
imshow3D(log(abs(squeeze(TWIST_Kspace(50,:,:,:,1)))))
title('K-space being acquired each measurement, Coil 1')

figure
imshow3D(log(abs(squeeze(unTWISTed_Kspace(50,:,:,:,1)))))
title('K-space being updated with each measurement, Coil 1')

figure
imshow3D(abs(squeeze(unTWISTed_IMspace(50,:,:,:,1))))

figure
imshow3D(log(abs(squeeze(Output_Coils_K(50,:,:,:,1)))))

imslice(abs(Output_Coils_IMspace))
title('Coil Images')

imslice(abs(SOS))
title('Sum of Squares')






