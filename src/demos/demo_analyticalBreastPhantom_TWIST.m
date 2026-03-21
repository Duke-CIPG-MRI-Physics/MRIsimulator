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


% Contrast parameters
rBW_HzPerPix = 570;
TR = (5.88E-3);
TE = 2.63E-3;

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*matrix_size_acquired(1);
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

%% Configure acquisition ordering and timing

pA = .04;
pB = .1;
shareMode = "forward";      % "forward" | "reverse" | "symmetric"
shareMethod = "single_anchor"; % use "dual_anchor" for symmetric sharing
shareTieBreaker = "future"; % symmetric tie-breaker: prefer later frame when equidistant
twistShareOptions = getTWISTShareOptions(struct( ...
    'shareMode', shareMode, ...
    'shareMethod', shareMethod, ...
    'shareTieBreaker', shareTieBreaker));

Num_Measurements = 5;
R = 1; %[2 3] %check if motion in breat sim is slowed down
PF_Factor = 1; %[6/8 6/8]

[Sampling_Table,TWIST_Timing] = Ultrafast_Sampling( ...
    matrix_size_acquired, FOV_acquired, pA, pB, Num_Measurements, TR, R, PF_Factor, ...
    twistShareOptions.mode, twistShareOptions.method);

%displaying region A
figure
regionA = getRegionA(matrix_size_acquired,FOV_acquired,pA,PF_Factor,R);
imshow(regionA)
title('Region A within Acquired Matrix')

%Construct Timing Information
TR_before_readout = TR-(matrix_size_acquired(1)*dt_s);
dwell_time_timepoints_within_TR = TR_before_readout+dt_s:dt_s:TR;

n_TRs_total = 1:(height(Sampling_Table)/matrix_size_acquired(1));

dwell_time_timepoints_absolute = dwell_time_timepoints_within_TR(:) + (n_TRs_total * TR);

Sampling_Table.Timing = dwell_time_timepoints_absolute(:);
clear n_TRs_total dwell_time_timepoints_absolute dwell_time_timepoints_within_TR

twistPlan = prepareTWISTViewSharingPlan( ...
    Sampling_Table, matrix_size_acquired, twistShareOptions.mode, twistShareOptions.tieBreaker);


%% Build 1D WORLD k-space grids
disp('Building WORLD k-space')
k_spatFreq_freq = computeKspaceGrid1D(FOV_acquired(1), matrix_size_acquired(1));
k_spatFreq_phase = computeKspaceGrid1D(FOV_acquired(2), matrix_size_acquired(2));
k_spatFreq_slice = computeKspaceGrid1D(FOV_acquired(3), matrix_size_acquired(3));
fps_to_xyz = zeros(1, 3);
fps_to_xyz(freq_phase_slice) = 1:3;

%% Construct Breast Phantom

% Update startInjectionTime_s to be relative to first frame ending
firstFrameMask = Sampling_Table.Frame == twistPlan.frameNumbers(1);
endOfFirstFrame = max(Sampling_Table.Timing(firstFrameMask));

% Injected contrast parameters
breastPhantomParams.startInjectionTime_s = breastPhantomParams.startInjectionTime_s + endOfFirstFrame;
breastPhantomParams.lesionArrivalDelay_s = 85;
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

%% Perform TWIST with streaming/bounded-buffer view sharing
maxChunkSize = 5000000;
nTimes = twistPlan.nFrames;
noiseSigma = 100;

padsize = matrix_size_complete(fps_to_xyz) - matrix_size_acquired(fps_to_xyz);
%% --- 8. Resolving Oversampling
crop_amount = matrix_size_complete-IMmatrix_crop_size;
margin = floor(crop_amount ./ 2);
cropRanges = { ...
    margin(1)+1 : margin(1)+IMmatrix_crop_size(1), ...
    margin(2)+1 : margin(2)+IMmatrix_crop_size(2), ...
    margin(3)+1 : margin(3)+IMmatrix_crop_size(3)};
voxel_volume = prod(voxel_size_mm);
final_IMspace = zeros([IMmatrix_crop_size, nTimes]);

switch twistShareOptions.mode
    case "forward"
        sharedKspace = [];
        for iTime = 1:nTimes
            fprintf('Reconstructing TWIST frame %d of %d (%.1f%% complete).\n', ...
                iTime, nTimes, (iTime - 1) / nTimes * 100);

            frameSamples = sampleAnalyticalTWISTFrame( ...
                phantom, Sampling_Table, twistPlan, iTime, ...
                k_spatFreq_freq, k_spatFreq_phase, k_spatFreq_slice, ...
                freq_phase_slice, noiseSigma, maxChunkSize);

            if isempty(sharedKspace)
                sharedKspace = zeros(matrix_size_acquired, 'like', frameSamples.aValues);
            end

            sharedKspace(twistPlan.aLinearIdx) = frameSamples.aValues;
            if twistPlan.frameIsFull(iTime)
                subsetsToUpdate = 1:twistPlan.nBSubsets;
            elseif twistPlan.acquiredSubsetByFrame(iTime) > 0
                subsetsToUpdate = twistPlan.acquiredSubsetByFrame(iTime);
            else
                subsetsToUpdate = [];
            end

            for iSubset = subsetsToUpdate
                sharedKspace(twistPlan.bLinearIdxBySubset{iSubset}) = frameSamples.bValues{iSubset};
            end

            final_IMspace(:,:,:,iTime) = reconstructTwistFrameImage( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, voxel_volume);
        end

    case "reverse"
        sharedKspace = [];
        for iTime = nTimes:-1:1
            fprintf('Reconstructing TWIST frame %d of %d (%.1f%% complete).\n', ...
                iTime, nTimes, (nTimes - iTime) / nTimes * 100);

            frameSamples = sampleAnalyticalTWISTFrame( ...
                phantom, Sampling_Table, twistPlan, iTime, ...
                k_spatFreq_freq, k_spatFreq_phase, k_spatFreq_slice, ...
                freq_phase_slice, noiseSigma, maxChunkSize);

            if isempty(sharedKspace)
                sharedKspace = zeros(matrix_size_acquired, 'like', frameSamples.aValues);
            end

            sharedKspace(twistPlan.aLinearIdx) = frameSamples.aValues;
            if twistPlan.frameIsFull(iTime)
                subsetsToUpdate = 1:twistPlan.nBSubsets;
            elseif twistPlan.acquiredSubsetByFrame(iTime) > 0
                subsetsToUpdate = twistPlan.acquiredSubsetByFrame(iTime);
            else
                subsetsToUpdate = [];
            end

            for iSubset = subsetsToUpdate
                sharedKspace(twistPlan.bLinearIdxBySubset{iSubset}) = frameSamples.bValues{iSubset};
            end

            final_IMspace(:,:,:,iTime) = reconstructTwistFrameImage( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, voxel_volume);
        end

    case "symmetric"
        frameBuffer = cell(1, nTimes);
        subsetReferenceCount = zeros(nTimes, twistPlan.nBSubsets);
        for iFrame = 1:nTimes
            for iSubset = 1:twistPlan.nBSubsets
                sourceFrame = twistPlan.sourceFrameMap(iFrame, iSubset);
                subsetReferenceCount(sourceFrame, iSubset) = subsetReferenceCount(sourceFrame, iSubset) + 1;
            end
        end

        nextFrameToReconstruct = 1;
        for iAcquire = 1:nTimes
            fprintf('Sampling TWIST frame %d of %d (%.1f%% complete).\n', ...
                iAcquire, nTimes, (iAcquire - 1) / nTimes * 100);

            frameBuffer{iAcquire} = sampleAnalyticalTWISTFrame( ...
                phantom, Sampling_Table, twistPlan, iAcquire, ...
                k_spatFreq_freq, k_spatFreq_phase, k_spatFreq_slice, ...
                freq_phase_slice, noiseSigma, maxChunkSize);

            while nextFrameToReconstruct <= nTimes && ...
                    twistPlan.readyFrameByFrame(nextFrameToReconstruct) <= iAcquire

                fprintf('Reconstructing TWIST frame %d of %d (%.1f%% complete).\n', ...
                    nextFrameToReconstruct, nTimes, ...
                    (nextFrameToReconstruct - 1) / nTimes * 100);

                currentSamples = frameBuffer{nextFrameToReconstruct};
                currentKspace = zeros(matrix_size_acquired, 'like', currentSamples.aValues);
                currentKspace(twistPlan.aLinearIdx) = currentSamples.aValues;

                for iSubset = 1:twistPlan.nBSubsets
                    sourceFrame = twistPlan.sourceFrameMap(nextFrameToReconstruct, iSubset);
                    currentKspace(twistPlan.bLinearIdxBySubset{iSubset}) = ...
                        frameBuffer{sourceFrame}.bValues{iSubset};

                    subsetReferenceCount(sourceFrame, iSubset) = ...
                        subsetReferenceCount(sourceFrame, iSubset) - 1;
                    if subsetReferenceCount(sourceFrame, iSubset) == 0
                        frameBuffer{sourceFrame}.bValues{iSubset} = [];
                    end
                end

                final_IMspace(:,:,:,nextFrameToReconstruct) = reconstructTwistFrameImage( ...
                    currentKspace, fps_to_xyz, padsize, cropRanges, voxel_volume);

                frameBuffer{nextFrameToReconstruct}.aValues = [];
                if frameSamplesEmpty(frameBuffer{nextFrameToReconstruct})
                    frameBuffer{nextFrameToReconstruct} = [];
                end

                for iSubset = 1:twistPlan.nBSubsets
                    sourceFrame = twistPlan.sourceFrameMap(nextFrameToReconstruct, iSubset);
                    if ~isempty(frameBuffer{sourceFrame}) && frameSamplesEmpty(frameBuffer{sourceFrame})
                        frameBuffer{sourceFrame} = [];
                    end
                end

                nextFrameToReconstruct = nextFrameToReconstruct + 1;
            end
        end

        if nextFrameToReconstruct <= nTimes
            error("demo_analyticalBreastPhantom_TWIST:IncompleteSymmetricReconstruction", ...
                "Symmetric TWIST reconstruction did not resolve all frame dependencies.");
        end
end


%% display phantom

phantom_magnitude = abs(final_IMspace);
imslice(squeeze(phantom_magnitude(:,:,:,:)));


% %% Contrast dynamics calculation
%
% %Ground truth
% figure;
% plot(Sampling_Table.Timing(1:1000:end),breastPhantomParams.lesionIntensityFunction(Sampling_Table.Timing(1:1000:end)) ...
%     + breastPhantomParams.breastIntensity);
%
% %convert TWIST frames to actual time, time for a whole frame is defined as
% %   moment when center of k-space is sampled.
%
% kspace_center = floor(matrix_size_acquired/2)+1;
%
% TWIST_frame_times = Sampling_Table.Timing((Sampling_Table.Frequency == kspace_center(1)) & ...
%        (Sampling_Table.("Row (phase)") == kspace_center(2)) & ...
%        (Sampling_Table.("Column (slice)") == kspace_center(3)));
%
%
% %TODO: build function to output lesion ROI
% lesion_center = [68,49,120]; %[freq,phase,slice] in final image
% lesion_radius = 6;
% [X, Y, Z] = ndgrid(1:IMmatrix_crop_size(1), 1:IMmatrix_crop_size(2), 1:IMmatrix_crop_size(3));
% squared_dist = (X - lesion_center(1)).^2 + (Y - lesion_center(2)).^2 + (Z - lesion_center(3)).^2;
% sphere_roi = squared_dist <= lesion_radius^2;
%
% data_reshaped = reshape(phantom_magnitude, [], size(phantom_magnitude,4));
% roi_flattened = sphere_roi(:);
% roi_data = data_reshaped(roi_flattened, :);
% contrast_values_measured = mean(roi_data, 1);
%
% hold on
% plot(TWIST_frame_times,abs(contrast_values_measured),'.-','MarkerSize',15)
% legend("Ground Truth","TWIST Measured")
% hold off
%
% title("Contrast Wash-in")
% xlabel("Time (s)")
% ylabel("Pixel Value")
%
% %%  Visualize ROI Overlay
% % 1. Select the slice to view (makes sense to use the lesion's Z-center)
% slice_to_view = lesion_center(3);
%
% % 2. Extract the final time frame for the background image
% final_time_idx = size(phantom_magnitude, 4);
% % Note: phantom_magnitude is 4D (freq, phase, slice, time)
% background_slice = phantom_magnitude(:, :, slice_to_view, final_time_idx);
%
% % 3. Extract the exact same slice from your 3D logical ROI mask
% roi_slice = sphere_roi(:, :, slice_to_view);
%
% % 4. Plotting
% figure;
% % Use imagesc for raw MRI data to automatically scale the display contrast
% imagesc(background_slice);
% colormap(gray);
% axis image; % Fixes the aspect ratio so the image isn't stretched
% axis off;   % Hides the axis ticks for a cleaner look
% title(sprintf('ROI Overlay on Slice %d (Final Time Frame)', slice_to_view));
%
% hold on;
% % Overlay the ROI as a red outline
% % The [0.5 0.5] tells contour to draw the line exactly at the logical boundary
% contour(roi_slice, [0.5 0.5], 'r', 'LineWidth', 2);
% hold off;
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

function croppedImage = reconstructTwistFrameImage(currentKspace, fps_to_xyz, padsize, cropRanges, voxel_volume)
% reconstructTwistFrameImage  Reconstruct and crop one TWIST image frame.

paddedKspace = padarray(permute(currentKspace, fps_to_xyz), 0.5 * padsize, 0);
currentImage = fftshift(ifftn(ifftshift(paddedKspace)));
currentImage = permute(currentImage, [2, 1, 3]);
currentImage = currentImage ./ voxel_volume;
croppedImage = currentImage(cropRanges{1}, cropRanges{2}, cropRanges{3});
end

function isEmpty = frameSamplesEmpty(frameSamples)
% frameSamplesEmpty  True when a buffered frame no longer stores raw data.

if isempty(frameSamples)
    isEmpty = true;
    return
end

if ~isempty(frameSamples.aValues)
    isEmpty = false;
    return
end

isEmpty = all(cellfun(@isempty, frameSamples.bValues));
end
