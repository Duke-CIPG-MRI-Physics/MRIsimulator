function [output] = Analytical_TWIST_Simulator(SimulationParameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% arguments (Input)
%     inputArg1
%     inputArg2
% end

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

[Sampling_Table,TWIST_Timing] = Ultrafast_Sampling(matrix_size_acquired,FOV_acquired,pA,pB,Num_Measurements,TR,R,PF_Factor);

k_idx_freq_pha_sli = [Sampling_Table.Frequency, Sampling_Table.("Row (phase)"), Sampling_Table.("Column (slice)")];

%Construct Timing Information
TR_before_readout = TR-(matrix_size_acquired(1)*dt_s); 
dwell_time_timepoints_within_TR = TR_before_readout+dt_s:dt_s:TR;

n_TRs_total = 1:(height(Sampling_Table)/matrix_size_acquired(1));

dwell_time_timepoints_absolute = dwell_time_timepoints_within_TR(:) + (n_TRs_total * TR);

Sampling_Table.Timing = dwell_time_timepoints_absolute(:);
clear n_TRs_total dwell_time_timepoints_absolute dwell_time_timepoints_within_TR


%% Build WORLD k-space grid and map to the TWIST ordering
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
breastPhantomParams.startInjectionTime_s = endOfFirstFrame;
breastPhantomParams.lesionIntensityFunction = @(t_s) calculateLesionEnhancement( ...
    t_s, breastPhantomParams, breastPhantomParams.lesionWashinType, ...
    breastPhantomParams.lesionWashoutType, breastPhantomParams.lesionKineticOverrides);


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
        
    % 1. Calculate current K-space Samples and put points in correct locations
    currentMask = (Sampling_Table.Frame == (iTime - 1));
    currentKspace = nan(matrix_size_acquired);
    
    currentIdx = sub2ind(matrix_size_acquired, ...
        Sampling_Table.Frequency(currentMask), ...
        Sampling_Table.("Row (phase)")(currentMask), ...
        Sampling_Table.("Column (slice)")(currentMask));
        
    currentKspace(currentIdx) = phantom.kspaceAtTime(...
        k_spatFreq_xyz(1, currentMask), ...
        k_spatFreq_xyz(2, currentMask), ...
        k_spatFreq_xyz(3, currentMask), ...
        Sampling_Table.Timing(currentMask)', ...
        maxChunkSize)';
    
    % 2. Apply View Sharing (Fill in missing k-space points if pB is active)
    if pB ~= 0 && iTime > 1
        missingKspaceData = isnan(currentKspace);
        currentKspace(missingKspaceData) = previousKspace(missingKspaceData);
    end
    
    % 3. Zero-fill unsampled periphery (Crucial for pB == 0)
    currentKspace(isnan(currentKspace)) = 0; 
    
    % 4. Reconstruct Image (Permute, zero-pad, shift, and IFFT)
    paddedKspace = padarray(permute(currentKspace, fps_to_xyz), 0.5 * padsize, 0);
    twistImage(:,:,:,iTime) = fftshift(ifftn(ifftshift(paddedKspace)));
    
    % 5. Prepare for the next TWIST frame
    if pB ~= 0
        previousKspace = currentKspace;
    end
end
clear k_spatFreq_xyz 

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

%% Contrast dynamics calculation


%convert TWIST frames to actual time, time for a whole frame is defined as
%   moment when center of k-space is sampled.

kspace_center = floor(matrix_size_acquired/2)+1;

TWIST_frame_times = (Sampling_Table.Timing((Sampling_Table.Frequency == kspace_center(1)) & ...
       (Sampling_Table.("Row (phase)") == kspace_center(2)) & ...
       (Sampling_Table.("Column (slice)") == kspace_center(3))))';


%TODO: build function to output lesion ROI
lesion_center = [68,49,120]; %[freq,phase,slice] in final image
lesion_radius = 6;
[X, Y, Z] = ndgrid(1:IMmatrix_crop_size(1), 1:IMmatrix_crop_size(2), 1:IMmatrix_crop_size(3));
squared_dist = (X - lesion_center(1)).^2 + (Y - lesion_center(2)).^2 + (Z - lesion_center(3)).^2;
sphere_roi = squared_dist <= lesion_radius^2;

data_reshaped = reshape(abs(final_IMspace), [], size(final_IMspace,4));
roi_flattened = sphere_roi(:);
roi_data = data_reshaped(roi_flattened, :);
contrast_values_measured = mean(roi_data, 1);


output.measured.contrast = contrast_values_measured;

% 1. Shift measured timepoints to be relative to lesion arrival
output.measured.timepoints = TWIST_frame_times - ...
    (breastPhantomParams.startInjectionTime_s+breastPhantomParams.lesionArrivalDelay_s);

% 2. Generate the simulated curve over a standardized relative time window
rel_time_vector = -20:300; % Adjust this window as needed for your wash-in/wash-out
abs_time_vector = rel_time_vector + ...
    (breastPhantomParams.startInjectionTime_s+breastPhantomParams.lesionArrivalDelay_s);

output.simulated.contrast = breastPhantomParams.lesionIntensityFunction(abs_time_vector) + breastPhantomParams.breastIntensity;
output.simulated.timepoints = rel_time_vector;
toc
end
