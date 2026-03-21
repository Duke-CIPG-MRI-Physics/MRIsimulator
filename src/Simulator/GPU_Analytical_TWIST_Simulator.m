function [output] = GPU_Analytical_TWIST_Simulator(SimulationParameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Unpacking inputs
tic
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
rBW_HzPerPix = SimulationParameters.MRIContrastParameters.rBW_HzPerPix;
TR = SimulationParameters.MRIContrastParameters.TR;
TE = SimulationParameters.MRIContrastParameters.TE;

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

phantom = BreastPhantom(breastPhantomParams);

%% ================= GPU ACCELERATION SETUP ================= 
k_spatFreq_freq_gpu = gpuArray(k_spatFreq_freq);
k_spatFreq_phase_gpu = gpuArray(k_spatFreq_phase);
k_spatFreq_slice_gpu = gpuArray(k_spatFreq_slice);
%% ==========================================================

%% Streaming TWIST reconstruction
maxChunkSize = 5000000;
nTimes = twistPlan.nFrames;
noiseSigma = 100;
padsize = matrix_size_complete(fps_to_xyz) - matrix_size_acquired(fps_to_xyz);

%% --- Resolving Oversampling on GPU
crop_amount = matrix_size_complete-IMmatrix_crop_size;
margin = floor(crop_amount ./ 2); 
cropRanges = { ...
    margin(1)+1 : margin(1)+IMmatrix_crop_size(1), ...
    margin(2)+1 : margin(2)+IMmatrix_crop_size(2), ...
    margin(3)+1 : margin(3)+IMmatrix_crop_size(3)};

%% Contrast dynamics calculation
%convert TWIST frames to actual time, time for a whole frame is defined as
%   moment when center of k-space is sampled.
TWIST_frame_times = twistPlan.frameCenterTimes_s;

%TODO: build function to output lesion ROI
lesion_center = [68,49,120]; %[freq,phase,slice] in final image
lesion_radius = 6;
[X, Y, Z] = ndgrid(1:IMmatrix_crop_size(1), 1:IMmatrix_crop_size(2), 1:IMmatrix_crop_size(3));

squared_dist = (X - lesion_center(1)).^2 + (Y - lesion_center(2)).^2 + (Z - lesion_center(3)).^2;
sphere_roi = squared_dist <= lesion_radius^2;
sphere_roi_gpu = gpuArray(sphere_roi);
voxel_volume = prod(voxel_size_mm);
contrast_values_measured = zeros(1, nTimes);

switch shareOptions.mode
    case "forward"
        sharedKspace = [];
        for iTime = 1:nTimes
            frameSamples = sampleAnalyticalTWISTFrame( ...
                phantom, Sampling_Table, twistPlan, iTime, ...
                k_spatFreq_freq_gpu, k_spatFreq_phase_gpu, k_spatFreq_slice_gpu, ...
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

            contrast_values_measured(iTime) = reconstructTwistRoiMeanGpu( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, sphere_roi_gpu, voxel_volume);
        end

    case "reverse"
        sharedKspace = [];
        for iTime = nTimes:-1:1
            frameSamples = sampleAnalyticalTWISTFrame( ...
                phantom, Sampling_Table, twistPlan, iTime, ...
                k_spatFreq_freq_gpu, k_spatFreq_phase_gpu, k_spatFreq_slice_gpu, ...
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

            contrast_values_measured(iTime) = reconstructTwistRoiMeanGpu( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, sphere_roi_gpu, voxel_volume);
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
            frameBuffer{iAcquire} = sampleAnalyticalTWISTFrame( ...
                phantom, Sampling_Table, twistPlan, iAcquire, ...
                k_spatFreq_freq_gpu, k_spatFreq_phase_gpu, k_spatFreq_slice_gpu, ...
                freq_phase_slice, noiseSigma, maxChunkSize);

            while nextFrameToReconstruct <= nTimes && ...
                    twistPlan.readyFrameByFrame(nextFrameToReconstruct) <= iAcquire

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

                contrast_values_measured(nextFrameToReconstruct) = reconstructTwistRoiMeanGpu( ...
                    currentKspace, fps_to_xyz, padsize, cropRanges, sphere_roi_gpu, voxel_volume);

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
            error("GPU_Analytical_TWIST_Simulator:IncompleteSymmetricReconstruction", ...
                "Symmetric TWIST reconstruction did not resolve all frame dependencies.");
        end
end

output.measured.contrast = contrast_values_measured;

% 1. Shift measured timepoints so t=0 is the injection start
output.measured.timepoints = TWIST_frame_times - breastPhantomParams.startInjectionTime_s;

% 2. Generate the simulated curve over a standardized relative time window
rel_time_vector = -20:300; % Adjust this window as needed for your wash-in/wash-out
abs_time_vector = rel_time_vector + breastPhantomParams.startInjectionTime_s;
output.simulated.contrast = breastPhantomParams.lesionIntensityFunction(abs_time_vector) + breastPhantomParams.breastIntensity;
output.simulated.timepoints = rel_time_vector;
toc
end

function roiMean = reconstructTwistRoiMeanGpu(currentKspace, fps_to_xyz, padsize, cropRanges, sphere_roi_gpu, voxel_volume)
% reconstructTwistRoiMeanGpu  Reconstruct one TWIST frame on the GPU and measure the ROI.

paddedKspace = padarray(permute(currentKspace, fps_to_xyz), 0.5 * padsize, 0);
currentImage = fftshift(ifftn(ifftshift(paddedKspace)));
currentImage = permute(currentImage, [2, 1, 3]);
currentImage = currentImage ./ voxel_volume;
currentImage = currentImage(cropRanges{1}, cropRanges{2}, cropRanges{3});

roiMean = gather(mean(abs(currentImage(sphere_roi_gpu))));
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
