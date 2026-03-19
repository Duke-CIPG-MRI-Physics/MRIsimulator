clear;
close all;
clc;

%% FOV and matrix size (scanner-style inputs)
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I
encodingFullStr = formatEncodingString(freq_phase_slice);
disp(encodingFullStr)

breastPhantomParams = createBreastPhantomParams();

%load('fast_scan_parameters.mat')
load('Breast_Ultrafast_scan_parameters.mat')
[FOV_acquired,matrix_size_complete,matrix_size_acquired,voxel_size_mm,nyquist_resolution_mm,IMmatrix_crop_size] =...
    convert_Siemens_parameters(scan_parameters);

cheat_factor = 3*2*(8/6)*(8/6);

% Contrast parameters
rBW_HzPerPix = 570*cheat_factor;
TR = (5.88E-3)/cheat_factor;
TE = 2.63E-3/cheat_factor;

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*matrix_size_acquired(1);
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

%% Configure acquisition ordering and timing

pA = .04;
pB = 0;

Num_Measurements = 10;
R = 1; %[2 3] 
PF_Factor = 1; %[6/8 6/8]

[Sampling_Table,TWIST_Timing] = Ultrafast_Sampling(matrix_size_acquired,FOV_acquired,pA,pB,Num_Measurements,TR,R,PF_Factor);

%displaying region A
figure
regionA = getRegionA(matrix_size_acquired,FOV_acquired,pA,PF_Factor,R);
imshow(regionA)
title('Region A within Acquired Matrix')

k_idx_freq_pha_sli = [Sampling_Table.Frequency, Sampling_Table.("Row (phase)"), Sampling_Table.("Column (slice)")];

%Construct Timing Information
TR_before_readout = TR-(matrix_size_acquired(1)*dt_s);
dwell_time_timepoints_within_TR = TR_before_readout+dt_s:dt_s:TR;

n_TRs_total = 1:(height(Sampling_Table)/matrix_size_acquired(1));

dwell_time_timepoints_absolute = dwell_time_timepoints_within_TR(:) + (n_TRs_total * TR);

Sampling_Table.Timing = dwell_time_timepoints_absolute(:);
clear n_TRs_total dwell_time_timepoints_absolute dwell_time_timepoints_within_TR


%% Build WORLD k-space grid and map to the TWIST ordering
disp('Building WORLD k-space')
k_spatFreq_freq = computeKspaceGrid1D(FOV_acquired(1), matrix_size_acquired(1));
k_spatFreq_phase = computeKspaceGrid1D(FOV_acquired(2), matrix_size_acquired(2));
k_spatFreq_slice = computeKspaceGrid1D(FOV_acquired(3), matrix_size_acquired(3));

k_spatFreq_freq_pha_sli = [k_spatFreq_freq(k_idx_freq_pha_sli(:, 1));
    k_spatFreq_phase(k_idx_freq_pha_sli(:, 2));
    k_spatFreq_slice(k_idx_freq_pha_sli(:, 3))];
[k_spatFreq_xyz, fps_to_xyz] = mapKspaceFpsToXyz(k_spatFreq_freq_pha_sli, freq_phase_slice);
clear k_spatFreq_freq_pha_sli k_spatFreq_freq k_spatFreq_phase k_spatFreq_slice k_idx_freq_pha_sli

%% Construct Breast Phantom

% Update startInjectionTime_s to be relative to first frame ending
endOfFirstFrame = max(Sampling_Table.Timing(Sampling_Table.Frame == 0));

% Injected contrast parameters
breastPhantomParams.startInjectionTime_s = breastPhantomParams.startInjectionTime_s + endOfFirstFrame;
breastPhantomParams.lesionArrivalDelay_s = 1;
breastPhantomParams.lesionWashinType = "instant";
breastPhantomParams.lesionWashoutType = "washout";
breastPhantomParams.lesionPeakEnhancement = 1.6;
breastPhantomParams.lesionBaselineDeltaIntensity = 0;
breastPhantomParams.lesionIntensityFunction = @(t_s) calculateLesionEnhancement( ...
    t_s, breastPhantomParams, breastPhantomParams.lesionWashinType, ...
    breastPhantomParams.lesionWashoutType, breastPhantomParams.lesionKineticOverrides);

sharedLesionIntensityFunction = breastPhantomParams.lesionIntensityFunction;
breastPhantomParams.lesions = [ ...
    struct('center_mm', [0, 0, 0], 'radius_mm', 10, 'intensityFunction', sharedLesionIntensityFunction), ...
    struct('center_mm', [28, 12, 16], 'radius_mm', 5, 'intensityFunction', sharedLesionIntensityFunction), ...
    struct('center_mm', [-16, -8, -14], 'radius_mm', 2.5, 'intensityFunction', sharedLesionIntensityFunction), ...
    struct('center_mm', [10, -17, 22], 'radius_mm', 1.25, 'intensityFunction', sharedLesionIntensityFunction)];

phantom = BreastPhantom(breastPhantomParams);

%% Perform TWIST, calculating two time frames at a time to minimize memory overhead
maxChunkSize = 5000000;
nTimes = Num_Measurements + 1;

% Preallocate the final image array
padsize = matrix_size_complete(fps_to_xyz) - matrix_size_acquired(fps_to_xyz);
twistImage = zeros([matrix_size_complete(fps_to_xyz), nTimes]);

% Initialize variable for view sharing
previousKspace = [];

for iTime = 1:nTimes
    fprintf('Reconstructing TWIST time %d of %d (%.1f%% complete).\n', ...
        iTime, nTimes, (iTime-1)/nTimes*100);

    % 1. Calculate current K-space Samples and put points in correct locations
    currentMask = (Sampling_Table.Frame == (iTime - 1));
    currentKspace = nan(matrix_size_acquired);

    currentIdx = sub2ind(matrix_size_acquired, ...
        Sampling_Table.Frequency(currentMask), ...
        Sampling_Table.("Row (phase)")(currentMask), ...
        Sampling_Table.("Column (slice)")(currentMask));

    temp_kspace = phantom.kspaceAtTime(...
        k_spatFreq_xyz(1, currentMask), ...
        k_spatFreq_xyz(2, currentMask), ...
        k_spatFreq_xyz(3, currentMask), ...
        Sampling_Table.Timing(currentMask)', ...
        maxChunkSize)';

    %Add noise
    sigma = 100;

    % Generate complex Gaussian noise
    % randn generates normal distribution N(0,1)
    noise = sigma * (randn(size(temp_kspace)) + 1i * randn(size(temp_kspace)));

    % Add to noiseless data
    temp_kspace = temp_kspace + noise;
    currentKspace(currentIdx) = temp_kspace;

    % 2. Apply View Sharing (Fill in missing k-space points)
    if  iTime > 1
        missingKspaceData = isnan(currentKspace);
        currentKspace(missingKspaceData) = previousKspace(missingKspaceData);
    end

    % 3. Zero-fill unsampled periphery (Crucial for pB == 0)
    currentKspace(isnan(currentKspace)) = 0;

    % 4. Reconstruct Image (Permute, zero-pad, shift, and IFFT)
    paddedKspace = padarray(permute(currentKspace, fps_to_xyz), 0.5 * padsize, 0);
    twistImage(:,:,:,iTime) = fftshift(ifftn(ifftshift(paddedKspace)));

    % 5. Prepare for the next TWIST frame
    previousKspace = currentKspace;

end
clear k_spatFreq_xyz temp_kspace

% Post-processing dimensions and intensity
twistImage = permute(twistImage, [2, 1, 3, 4]);

voxel_volume = prod(voxel_size_mm);
twistImage = twistImage ./ voxel_volume;
%% --- 8. Resolving Oversampling
crop_amount = matrix_size_complete-IMmatrix_crop_size;
margin = floor(crop_amount ./ 2);

final_IMspace = twistImage(...
    margin(1)+1 : margin(1)+IMmatrix_crop_size(1), ...
    margin(2)+1 : margin(2)+IMmatrix_crop_size(2), ...
    margin(3)+1 : margin(3)+IMmatrix_crop_size(3),...
    :);


%% display phantom

phantom_magnitude = abs(final_IMspace);
imslice(squeeze(phantom_magnitude(:,:,:,:)));


%% Contrast dynamics calculation

%Ground truth
figure;
plot(Sampling_Table.Timing(1:1000:end),breastPhantomParams.lesionIntensityFunction(Sampling_Table.Timing(1:1000:end)) ...
    + breastPhantomParams.breastIntensity);

%convert TWIST frames to actual time, time for a whole frame is defined as
%   moment when center of k-space is sampled.

kspace_center = floor(matrix_size_acquired/2)+1;

TWIST_frame_times = Sampling_Table.Timing((Sampling_Table.Frequency == kspace_center(1)) & ...
       (Sampling_Table.("Row (phase)") == kspace_center(2)) & ...
       (Sampling_Table.("Column (slice)") == kspace_center(3)));


%TODO: build function to output lesion ROI
lesion_center = [68,49,120]; %[freq,phase,slice] in final image
lesion_radius = 6;
[X, Y, Z] = ndgrid(1:IMmatrix_crop_size(1), 1:IMmatrix_crop_size(2), 1:IMmatrix_crop_size(3));
squared_dist = (X - lesion_center(1)).^2 + (Y - lesion_center(2)).^2 + (Z - lesion_center(3)).^2;
sphere_roi = squared_dist <= lesion_radius^2;

data_reshaped = reshape(phantom_magnitude, [], size(phantom_magnitude,4));
roi_flattened = sphere_roi(:);
roi_data = data_reshaped(roi_flattened, :);
contrast_values_measured = mean(roi_data, 1);

hold on
plot(TWIST_frame_times,abs(contrast_values_measured),'.-','MarkerSize',15)
legend("Ground Truth","TWIST Measured")
hold off

title("Contrast Wash-in")
xlabel("Time (s)")
ylabel("Pixel Value")

%%  Visualize ROI Overlay
% 1. Select the slice to view (makes sense to use the lesion's Z-center)
slice_to_view = lesion_center(3);

% 2. Extract the final time frame for the background image
final_time_idx = size(phantom_magnitude, 4);
% Note: phantom_magnitude is 4D (freq, phase, slice, time)
background_slice = phantom_magnitude(:, :, slice_to_view, final_time_idx);

% 3. Extract the exact same slice from your 3D logical ROI mask
roi_slice = sphere_roi(:, :, slice_to_view);

% 4. Plotting
figure;
% Use imagesc for raw MRI data to automatically scale the display contrast
imagesc(background_slice);
colormap(gray);
axis image; % Fixes the aspect ratio so the image isn't stretched
axis off;   % Hides the axis ticks for a cleaner look
title(sprintf('ROI Overlay on Slice %d (Final Time Frame)', slice_to_view));

hold on;
% Overlay the ROI as a red outline
% The [0.5 0.5] tells contour to draw the line exactly at the logical boundary
contour(roi_slice, [0.5 0.5], 'r', 'LineWidth', 2);
hold off;
%
% %% Saving output
% save_ask = input('Save output?: (y/n)','s');
%
% if strcmpi(save_ask, 'y')
%
%     %measuring phantom size
%     output_bytes = whos("phantom_magnitude");
%     output_bytes = output_bytes.bytes;
%
%     %creating structure with all phantom information
%
%     % Outputs
%     phantom_simulated.outputs.phantom_magnitude = phantom_magnitude;
%     phantom_simulated.outputs.timing            = TWIST_frame_times;
%
%     % Inputs (Consolidating everything under .inputs)
%     phantom_simulated.inputs.breastPhantomParams = breastPhantomParams;
%     phantom_simulated.inputs.scan_parameters     = scan_parameters;
%     phantom_simulated.inputs.TR                  = TR;
%     phantom_simulated.inputs.TE                  = TE;
%     phantom_simulated.inputs.rBW_HzPerPix        = rBW_HzPerPix;
%
%     % Nested TWIST parameters
%     phantom_simulated.inputs.TWIST.pA            = pA;
%     phantom_simulated.inputs.TWIST.pB            = 1/Nb;
%     phantom_simulated.inputs.TWIST.Time_Measured = Time_Measured;
%
%     % Nested Undersampling parameters
%     phantom_simulated.inputs.Undersampling.GRAPPA_R  = R;
%     phantom_simulated.inputs.Undersampling.PF_Factor = PF_Factor;
%
%
%     fprintf('Input desired filename, file will be saved as <filename>.mat\n')
%     filename = input(':','s');
%
%     if output_bytes >= 1.99e9
%         fprintf('Saving using v7.3...\n')
%         save(filename,'phantom_simulated','-v7.3')
%         fprintf('File saved as %s.mat\n, ', filename);
%
%     else
%         fprintf('Saving using v7...\n')
%         save(filename,'phantom_simulated','-v7')
%         fprintf('File saved as %s.mat\n', filename);
%
%     end
%
% else
%     fprintf('Output not saved')
%
% end
