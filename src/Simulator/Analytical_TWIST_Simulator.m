function [output] = Analytical_TWIST_Simulator(SimulationParameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% arguments (Input)
%     inputArg1
%     inputArg2
% end

%% Unpacking inputs
tic

cheat_factor = 2*3*(8/6)*(8/6);

%Scan Parameters
scan_parameters = SimulationParameters.ScanParameters;

% Lesion parameters
breastPhantomParams = createBreastPhantomParams();
breastPhantomParams.lesionArrivalDelay_s = SimulationParameters.LesionParameters.lesionArrivalDelay_s;
breastPhantomParams.lesionWashinType = SimulationParameters.LesionParameters.lesionWashinType;
breastPhantomParams.lesionWashoutType = SimulationParameters.LesionParameters.lesionWashoutType;
breastPhantomParams.lesionPeakEnhancement = SimulationParameters.LesionParameters.lesionPeakEnhancement;
breastPhantomParams.lesionBaselineDeltaIntensity = SimulationParameters.LesionParameters.lesionBaselineDeltaIntensity;


% MRI contrast parameters
rBW_HzPerPix = SimulationParameters.MRIContrastParameters.rBW_HzPerPix*cheat_factor;
TR = SimulationParameters.MRIContrastParameters.TR/cheat_factor;
TE = SimulationParameters.MRIContrastParameters.TE/cheat_factor;

% Ultrafast parameters
pA = SimulationParameters.TWIST.pA;
pB = SimulationParameters.TWIST.pB;
Num_Measurements = SimulationParameters.TWIST.N_measurements;
shareOptions = getTWISTShareOptions(SimulationParameters.TWIST);

R = SimulationParameters.ParallelImaging.GRAPPA_R;
PF_Factor = SimulationParameters.ParallelImaging.PF_Factor;


%% FOV and matrix size (scanner-style inputs)
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I 
encodingFullStr = formatEncodingString(freq_phase_slice);


[FOV_acquired,matrix_size_complete,matrix_size_acquired,voxel_size_mm,nyquist_resolution_mm,IMmatrix_crop_size] =...
    convert_Siemens_parameters(scan_parameters);

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*matrix_size_acquired(1);  
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

%% Configure acquisition ordering and timing

[Sampling_Table,TWIST_Timing] = Ultrafast_Sampling( ...
    matrix_size_acquired, FOV_acquired, pA, pB, Num_Measurements, TR, R, PF_Factor, ...
    shareOptions.mode, shareOptions.method);

%Construct Timing Information
TR_before_readout = TR-(matrix_size_acquired(1)*dt_s); 
dwell_time_timepoints_within_TR = TR_before_readout+dt_s:dt_s:TR;

n_TRs_total = 1:(height(Sampling_Table)/matrix_size_acquired(1));

dwell_time_timepoints_absolute = dwell_time_timepoints_within_TR(:) + (n_TRs_total * TR);

Sampling_Table.Timing = dwell_time_timepoints_absolute(:);
clear n_TRs_total dwell_time_timepoints_absolute dwell_time_timepoints_within_TR

twistPlan = prepareTWISTViewSharingPlan( ...
    Sampling_Table, matrix_size_acquired, shareOptions.mode, shareOptions.tieBreaker);


%% Build 1D WORLD k-space grids
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
breastPhantomParams.lesionIntensityFunction = @(t_s) calculateLesionEnhancement( ...
    t_s, breastPhantomParams, breastPhantomParams.lesionWashinType, ...
    breastPhantomParams.lesionWashoutType, breastPhantomParams.lesionKineticOverrides);

sharedLesionIntensityFunction = breastPhantomParams.lesionIntensityFunction;
breastPhantomParams.lesions = [ ...
    struct('center_mm', [-130*voxel_size_mm(1), 0, 0], 'radius_mm', 8, 'intensityFunction', sharedLesionIntensityFunction), ...
    struct('center_mm', [18*voxel_size_mm(1), 30*voxel_size_mm(2), 16*voxel_size_mm(3)], 'radius_mm', 4, 'intensityFunction', sharedLesionIntensityFunction), ...
    struct('center_mm', [-10*voxel_size_mm(1), -10*voxel_size_mm(2), -14*voxel_size_mm(3)], 'radius_mm', 2, 'intensityFunction', sharedLesionIntensityFunction), ...
    struct('center_mm', [6*voxel_size_mm(1), -20*voxel_size_mm(2), 22]*voxel_size_mm(3), 'radius_mm', 1, 'intensityFunction', sharedLesionIntensityFunction)];

phantom = BreastPhantom(breastPhantomParams);

%% Streaming TWIST reconstruction
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
    margin(3)+1 : margin(3)+IMmatrix_crop_size(3),...
    :);

final_IMspace_mag = abs(final_IMspace);
%% Contrast dynamics calculation


%convert TWIST frames to actual time, time for a whole frame is defined as
%   moment when center of k-space is sampled.
TWIST_frame_times = twistPlan.frameCenterTimes_s;

%TODO: build function to output lesion ROI
lesion_centers = [68,179,121; 38,31,105; 78,59,135; 88,43,99]; %[freq,phase,slice] in final image
lesion_radii = [5; 2.5; 1; 0];
num_lesions = length(lesion_radii);
num_volumes = size(final_IMspace_mag, 4); % Assuming 4D data (e.g., multiple echoes/timepoints)

% 2. Create coordinate grid
[X, Y, Z] = ndgrid(1:IMmatrix_crop_size(1), 1:IMmatrix_crop_size(2), 1:IMmatrix_crop_size(3));
squared_dist = (X - lesion_center(1)).^2 + (Y - lesion_center(2)).^2 + (Z - lesion_center(3)).^2;
sphere_roi = squared_dist <= lesion_radius^2;
voxel_volume = prod(voxel_size_mm);
contrast_values_measured = zeros(1, nTimes);

switch shareOptions.mode
    case "forward"
        sharedKspace = [];
        for iTime = 1:nTimes
            fprintf('Reconstructing TWIST frame %d of %d (%.1f%% complete).\n', ...
                iTime, nTimes, (iTime - 1) / nTimes * 100);

% 3. Reshape the source data once: rows = voxels, columns = volumes
data_reshaped = reshape(final_IMspace_mag, [], num_volumes);

% Pre-allocate outputs for speed and memory efficiency
lesion_roi = false(IMmatrix_crop_size(1), IMmatrix_crop_size(2), IMmatrix_crop_size(3), num_lesions);
contrast_values_measured = zeros(num_lesions, num_volumes);

% 4. Loop through each lesion
for ii = 1:num_lesions
    % Calculate distance from this specific center
    squared_dists = (X - lesion_centers(ii,1)).^2 + ...
                    (Y - lesion_centers(ii,2)).^2 + ...
                    (Z - lesion_centers(ii,3)).^2;
    
    % Create the 3D mask for this specific lesion
    current_mask = squared_dists <= lesion_radii(ii)^2;
    
    % Save it into your 4D stack in case you need to visualize it later
    lesion_roi(:,:,:,ii) = current_mask;
    
    % Flatten the mask to 1D to match data_reshaped
    roi_flattened = current_mask(:);
    
    % Extract the data just for this lesion mask
    roi_data = data_reshaped(roi_flattened, :);
    
    % Calculate the mean across the ROI (dimension 1) and store it
    contrast_values_measured(ii, :) = mean(roi_data, 1);

end

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

            contrast_values_measured(iTime) = reconstructTwistRoiMean( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, sphere_roi, voxel_volume);
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

            contrast_values_measured(iTime) = reconstructTwistRoiMean( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, sphere_roi, voxel_volume);
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

                contrast_values_measured(nextFrameToReconstruct) = reconstructTwistRoiMean( ...
                    currentKspace, fps_to_xyz, padsize, cropRanges, sphere_roi, voxel_volume);

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
            error("Analytical_TWIST_Simulator:IncompleteSymmetricReconstruction", ...
                "Symmetric TWIST reconstruction did not resolve all frame dependencies.");
        end
end

output.measured.contrast_L = contrast_values_measured(1,:);
output.measured.contrast_M = contrast_values_measured(2,:);
output.measured.contrast_S = contrast_values_measured(3,:);
output.measured.contrast_XS = contrast_values_measured(4,:);

% 1. Shift measured timepoints to be relative to lesion arrival
output.measured.timepoints = TWIST_frame_times - ...
    (breastPhantomParams.startInjectionTime_s+breastPhantomParams.lesionArrivalDelay_s);

% 2. Generate the simulated curve over a standardized relative time window
rel_time_vector = -2:.1:60; % Adjust this window as needed for your wash-in/wash-out
abs_time_vector = rel_time_vector + ...
    (breastPhantomParams.startInjectionTime_s+breastPhantomParams.lesionArrivalDelay_s);

output.simulated.contrast = breastPhantomParams.lesionIntensityFunction(abs_time_vector) + breastPhantomParams.breastIntensity;
output.simulated.timepoints = rel_time_vector;
toc
end

function roiMean = reconstructTwistRoiMean(currentKspace, fps_to_xyz, padsize, cropRanges, sphere_roi, voxel_volume)
% reconstructTwistRoiMean  Reconstruct one TWIST frame and measure the ROI.

paddedKspace = padarray(permute(currentKspace, fps_to_xyz), 0.5 * padsize, 0);
currentImage = fftshift(ifftn(ifftshift(paddedKspace)));
currentImage = permute(currentImage, [2, 1, 3]);
currentImage = currentImage ./ voxel_volume;
currentImage = currentImage(cropRanges{1}, cropRanges{2}, cropRanges{3});

roiValues = abs(currentImage(sphere_roi));
roiMean = mean(roiValues(:));
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
