function output = Analytical_TWIST_Simulator(SimulationParameters)
% Analytical_TWIST_Simulator  Simulate analytical TWIST sampling of the breast phantom.
%   output = Analytical_TWIST_Simulator(SimulationParameters) evaluates the
%   analytical breast phantom directly in k-space using the configured TWIST
%   ordering and view-sharing mode, then reports lesion ROI kinetics.
%
%   Inputs
%   ------
%   SimulationParameters : struct
%       Simulation configuration created by Make_Simulation_Parameters.m.
%
%   Output
%   ------
%   output : struct
%       output.measured contains ROI curves for the reconstructed TWIST
%       frames and output.simulated contains the lesion ground-truth curve.

%% Unpack inputs
tic

cheat_factor = 2 * 3 * (8 / 6) * (8 / 6);

scan_parameters = SimulationParameters.ScanParameters;

breastPhantomParams = createBreastPhantomParams();
breastPhantomParams.lesionArrivalDelay_s = ...
    SimulationParameters.LesionParameters.lesionArrivalDelay_s;
breastPhantomParams.lesionWashinType = ...
    SimulationParameters.LesionParameters.lesionWashinType;
breastPhantomParams.lesionWashoutType = ...
    SimulationParameters.LesionParameters.lesionWashoutType;
breastPhantomParams.lesionPeakEnhancement = ...
    SimulationParameters.LesionParameters.lesionPeakEnhancement;
breastPhantomParams.lesionBaselineDeltaIntensity = ...
    SimulationParameters.LesionParameters.lesionBaselineDeltaIntensity;

rBW_HzPerPix = SimulationParameters.MRIContrastParameters.rBW_HzPerPix * cheat_factor;
TR = SimulationParameters.MRIContrastParameters.TR / cheat_factor;
TE = SimulationParameters.MRIContrastParameters.TE / cheat_factor; %#ok<NASGU>

pA = SimulationParameters.TWIST.pA;
pB = SimulationParameters.TWIST.pB;
Num_Measurements = SimulationParameters.TWIST.N_measurements;
shareOptions = getTWISTShareOptions(SimulationParameters.TWIST);
orderingOptions = getTWISTOrderingOptions(SimulationParameters.TWIST);

R = SimulationParameters.ParallelImaging.GRAPPA_R;
PF_Factor = SimulationParameters.ParallelImaging.PF_Factor;
if isscalar(R)
    R = [R, R];
end
if isscalar(PF_Factor)
    PF_Factor = [PF_Factor, PF_Factor];
end

%% FOV and matrix size
freq_phase_slice = [2 1 3]; % 1 = R/L, 2 = A/P, 3 = S/I

[FOV_acquired, matrix_size_complete, matrix_size_acquired, voxel_size_mm, ~, IMmatrix_crop_size] = ...
    convert_Siemens_parameters(scan_parameters);

rBW_Hz = rBW_HzPerPix * matrix_size_acquired(1);
dt_s = 1 / rBW_Hz;

%% Configure acquisition ordering and timing
[Sampling_Table, ~] = Ultrafast_Sampling( ...
    matrix_size_acquired, FOV_acquired, pA, pB, Num_Measurements, TR, R, PF_Factor, ...
    shareOptions.mode, shareOptions.method, orderingOptions);

TR_before_readout = TR - (matrix_size_acquired(1) * dt_s);
dwell_time_timepoints_within_TR = TR_before_readout + dt_s:dt_s:TR;
n_TRs_total = 1:(height(Sampling_Table) / matrix_size_acquired(1));
dwell_time_timepoints_absolute = dwell_time_timepoints_within_TR(:) + (n_TRs_total * TR);

Sampling_Table.Timing = dwell_time_timepoints_absolute(:);

twistPlan = prepareTWISTViewSharingPlan( ...
    Sampling_Table, matrix_size_acquired, shareOptions.mode, shareOptions.tieBreaker);

%% Build 1D WORLD k-space grids
k_spatFreq_freq = computeKspaceGrid1D(FOV_acquired(1), matrix_size_acquired(1));
k_spatFreq_phase = computeKspaceGrid1D(FOV_acquired(2), matrix_size_acquired(2));
k_spatFreq_slice = computeKspaceGrid1D(FOV_acquired(3), matrix_size_acquired(3));
fps_to_xyz = zeros(1, 3);
fps_to_xyz(freq_phase_slice) = 1:3;

%% Construct breast phantom
firstFrameMask = Sampling_Table.Frame == twistPlan.frameNumbers(1);
endOfFirstFrame = max(Sampling_Table.Timing(firstFrameMask));

breastPhantomParams.startInjectionTime_s = ...
    breastPhantomParams.startInjectionTime_s + endOfFirstFrame;
breastPhantomParams.lesionIntensityFunction = @(t_s) calculateLesionEnhancement( ...
    t_s, breastPhantomParams, breastPhantomParams.lesionWashinType, ...
    breastPhantomParams.lesionWashoutType, breastPhantomParams.lesionKineticOverrides);

sharedLesionIntensityFunction = breastPhantomParams.lesionIntensityFunction;
breastPhantomParams.lesions = createTwistLesionDefinitions(sharedLesionIntensityFunction);

phantom = BreastPhantom(breastPhantomParams);

%% Precompute reconstruction metadata
maxChunkSize = 5000000;
nTimes = twistPlan.nFrames;
noiseSigma = 100;
padsize = matrix_size_complete(fps_to_xyz) - matrix_size_acquired(fps_to_xyz);

crop_amount = matrix_size_complete - IMmatrix_crop_size;
margin = floor(crop_amount ./ 2);
cropRanges = { ...
    margin(1) + 1 : margin(1) + IMmatrix_crop_size(1), ...
    margin(2) + 1 : margin(2) + IMmatrix_crop_size(2), ...
    margin(3) + 1 : margin(3) + IMmatrix_crop_size(3)};

voxel_volume = prod(voxel_size_mm);
TWIST_frame_times = twistPlan.frameCenterTimes_s;
lesionBorder_mm = 1;
roiGridSpec = struct( ...
    'FOV_acquired_mm', FOV_acquired, ...
    'matrix_size_complete', matrix_size_complete, ...
    'IMmatrix_crop_size', IMmatrix_crop_size, ...
    'freq_phase_slice', freq_phase_slice);
[~, roiInfo] = phantom.buildLesionRoiMasks(roiGridSpec, lesionBorder_mm);
roiLinearIdxByLesion = roiInfo.linearIdxByLesion;
contrast_values_measured = zeros(numel(roiLinearIdxByLesion), nTimes);

%% Streaming TWIST reconstruction
switch shareOptions.mode
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
            subsetsToUpdate = getFrameSubsetUpdates(twistPlan, iTime);
            for iSubset = subsetsToUpdate
                sharedKspace(twistPlan.bLinearIdxBySubset{iSubset}) = frameSamples.bValues{iSubset};
            end

            contrast_values_measured(:, iTime) = reconstructTwistRoiMeans( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, roiLinearIdxByLesion, voxel_volume);
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
            subsetsToUpdate = getFrameSubsetUpdates(twistPlan, iTime);
            for iSubset = subsetsToUpdate
                sharedKspace(twistPlan.bLinearIdxBySubset{iSubset}) = frameSamples.bValues{iSubset};
            end

            contrast_values_measured(:, iTime) = reconstructTwistRoiMeans( ...
                sharedKspace, fps_to_xyz, padsize, cropRanges, roiLinearIdxByLesion, voxel_volume);
        end

    case "symmetric"
        frameBuffer = cell(1, nTimes);
        subsetReferenceCount = zeros(nTimes, twistPlan.nBSubsets);

        for iFrame = 1:nTimes
            for iSubset = 1:twistPlan.nBSubsets
                sourceFrame = twistPlan.sourceFrameMap(iFrame, iSubset);
                subsetReferenceCount(sourceFrame, iSubset) = ...
                    subsetReferenceCount(sourceFrame, iSubset) + 1;
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

                contrast_values_measured(:, nextFrameToReconstruct) = reconstructTwistRoiMeans( ...
                    currentKspace, fps_to_xyz, padsize, cropRanges, roiLinearIdxByLesion, voxel_volume);

                frameBuffer{nextFrameToReconstruct}.aValues = [];
                if frameSamplesEmpty(frameBuffer{nextFrameToReconstruct})
                    frameBuffer{nextFrameToReconstruct} = [];
                end

                for iSubset = 1:twistPlan.nBSubsets
                    sourceFrame = twistPlan.sourceFrameMap(nextFrameToReconstruct, iSubset);
                    if ~isempty(frameBuffer{sourceFrame}) && ...
                            frameSamplesEmpty(frameBuffer{sourceFrame})
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

%% Assemble outputs
output.measured.contrast_L = contrast_values_measured(1, :);
output.measured.contrast_M = contrast_values_measured(2, :);
output.measured.contrast_S = contrast_values_measured(3, :);
output.measured.contrast_XS = contrast_values_measured(4, :);
output.measured.timepoints = TWIST_frame_times - ...
    (breastPhantomParams.startInjectionTime_s + breastPhantomParams.lesionArrivalDelay_s);

rel_time_vector = -2:0.1:60;
abs_time_vector = rel_time_vector + ...
    (breastPhantomParams.startInjectionTime_s + breastPhantomParams.lesionArrivalDelay_s);

output.simulated.contrast = ...
    breastPhantomParams.lesionIntensityFunction(abs_time_vector) + ...
    breastPhantomParams.breastIntensity;
output.simulated.timepoints = rel_time_vector;
toc
end

function subsetsToUpdate = getFrameSubsetUpdates(twistPlan, frameIndex)
% getFrameSubsetUpdates  Return the Bj subsets acquired in one frame.

if twistPlan.frameIsFull(frameIndex)
    subsetsToUpdate = 1:twistPlan.nBSubsets;
elseif twistPlan.acquiredSubsetByFrame(frameIndex) > 0
    subsetsToUpdate = twistPlan.acquiredSubsetByFrame(frameIndex);
else
    subsetsToUpdate = [];
end
end

function roiMeans = reconstructTwistRoiMeans( ...
    currentKspace, fps_to_xyz, padsize, cropRanges, roiLinearIdxByLesion, voxel_volume)
% reconstructTwistRoiMeans  Reconstruct one TWIST frame and measure lesion ROIs.
%   Uses the raw [freq phase slice] storage order returned by the TWIST IFFT.

paddedKspace = padarray(permute(currentKspace, fps_to_xyz), 0.5 * padsize, 0);
currentImage = fftshift(ifftn(ifftshift(paddedKspace)));
currentImage = permute(currentImage, [2, 1, 3]);
currentImage = currentImage ./ voxel_volume;
currentImage = currentImage(cropRanges{1}, cropRanges{2}, cropRanges{3});
currentMagnitude = abs(currentImage);

nLesions = numel(roiLinearIdxByLesion);
roiMeans = zeros(nLesions, 1);
for iLesion = 1:nLesions
    roiValues = currentMagnitude(roiLinearIdxByLesion{iLesion});
    roiMeans(iLesion) = mean(roiValues(:));
end
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
