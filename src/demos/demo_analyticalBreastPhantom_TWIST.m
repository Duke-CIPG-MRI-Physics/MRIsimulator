clear all; 
close all; 
clc; 

%% 1) FOV and matrix size (scanner-style inputs)
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I 
encodingFullStr = formatEncodingString(freq_phase_slice);
disp(encodingFullStr)



load('fast_scan_parameters.mat') 
%load('base_scan_parameters.mat')
[FOV_acquired,matrix_size_complete,matrix_size_acquired,voxel_size_mm,nyquist_resolution_mm,IMmatrix_crop_size] = convert_Siemens_parameters(scan_parameters);

% Contrast parameters
rBW_HzPerPix = 570;
TR = (5.88E-3);  
TE = 2.63E-3;

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*matrix_size_acquired(1);  
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

%% 2) Configure acquisition ordering and timing

pA = 0.05;
Nb = 10;
Time_Measured = 300; %sec
R = 1; %[2 3]
PF_Factor = 1; %[6/8 6/8]

Sampling_Table = Ultrafast_Sampling(matrix_size_acquired,pA,Nb,Time_Measured,TR,R,PF_Factor);

k_idx_fps = [Sampling_Table.Frequency, Sampling_Table.("Row (phase)"), Sampling_Table.("Column (slice)")];

%% 3) Build WORLD k-space grid and map to the TWIST ordering
disp('Building WORLD k-space')
k_spatFreq_freq = computeKspaceGrid1D(FOV_acquired(1), matrix_size_acquired(1));
k_spatFreq_phase = computeKspaceGrid1D(FOV_acquired(2), matrix_size_acquired(2));
k_spatFreq_slice = computeKspaceGrid1D(FOV_acquired(3), matrix_size_acquired(3));

k_spatFreq_fps = [k_spatFreq_freq(k_idx_fps(:, 1)); 
    k_spatFreq_phase(k_idx_fps(:, 2)); 
    k_spatFreq_slice(k_idx_fps(:, 3))];
[k_spatFreq_xyz, fps_to_xyz] = mapKspaceFpsToXyz(k_spatFreq_fps, freq_phase_slice);
clear k_spatFreq_fps k_spatFreq_freq k_spatFreq_phase k_spatFreq_slice k_idx_fps

%% 5) Construct the breast phantom with the embedded enhancing vessel
t_PE = TR-(matrix_size_acquired(1)*dt_s); %TODO: rename these variables for clarity
t_s = t_PE+dt_s:dt_s:TR;
TR_counts = 1:height(Sampling_Table)/matrix_size_acquired(1);
matrix_result = t_s(:) + (TR_counts * TR);
Sampling_Table.Timing = matrix_result(:);
clear TR_counts t_s matrix_result

% % Force Timing to be increasing
% [sortedT,sortIdx] = sort(Sampling_Table.Timing);
% phantom = BreastPhantom(sortedT);

breastPhantomParams = createBreastPhantomParams();

%% Contrast timing visualization
timing_s = Sampling_Table.Timing;
contrastLength_mm = calculatePlugFlowInVessels(timing_s, breastPhantomParams);

figure('Name', 'TWIST Contrast Plug Flow and Frame Timing');
yyaxis left
plot(timing_s, contrastLength_mm, 'LineWidth', 1.5);
ylabel('Contrast length [mm]')
xlabel('Time [s]')

yyaxis right
frameNumbers = unique(Sampling_Table.Bj);
frameNumbers = frameNumbers(frameNumbers > 0);
for frameIdx = 1:numel(frameNumbers)
    frameNumber = frameNumbers(frameIdx);
    frameTimes_s = timing_s(Sampling_Table.Bj == frameNumber);
    if isempty(frameTimes_s)
        continue;
    end
    frameStart_s = min(frameTimes_s);
    frameEnd_s = max(frameTimes_s);
    rectTime_s = [frameStart_s frameStart_s frameEnd_s frameEnd_s];
    rectAmp = [0 frameNumber frameNumber 0];
    plot(rectTime_s, rectAmp, 'LineWidth', 1.0);
    hold on
end
ylabel('TWIST frame index')
ylim([0 max(frameNumbers) + 1])
grid on

phantom = BreastPhantom(breastPhantomParams);

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');
kspace = nan([matrix_size_acquired,max(Sampling_Table.Bj)+1]);

% Sample phantom at measured time points
kspace(Sampling_Table.("Linear Index")) = ...
    phantom.kspaceAtTime(k_spatFreq_xyz(1, :), k_spatFreq_xyz(2, :), ...
    k_spatFreq_xyz(3, :),Sampling_Table.Timing')';

clear k_spatFreq_xyz Sampling_Table phantom

%% ---  6. Updating K-Space from each measurement to undo TWIST
fprintf('\nUndoing TWIST...')

% Loop only through timepoints (vectorized over coils)
for ii_timepoint = 2:size(kspace,4)
    % % 1. Start with the data measured for the current time point
    % currentData = kspace(:,:,:,ii_timepoint,:);
    % 
    % % 2. We will use the previous data to fill in unmeasured data
    % previousData = kspace(:,:,:,ii_timepoint-1,:);
    % 
    % % 3. Any nan values were not calculated this time frame, so fill them
    % % in from the previous time frame
    % currentData(isnan(currentData)) = previousData(isnan(currentData));
    % 
    % % 4. Assign the updated slice back into the main array.
    % kspace(:,:,:,ii_timepoint,:) = currentData;

    % More memory efficient calculation of the steps above
    kspace(:,:,:,ii_timepoint,:) = kspace(:,:,:,ii_timepoint,:).*(~isnan(kspace(:,:,:,ii_timepoint,:))) ...
        + kspace(:,:,:,ii_timepoint-1,:).*(~isnan(kspace(:,:,:,ii_timepoint,:)));
end

clear currentData previousData 

%% FFT and zero padd, one 3D volume at a time to reduce memory spikes
kspace = permute(kspace, [fps_to_xyz, 4 5]);
twistImage = zeros([matrix_size_complete(fps_to_xyz) size(kspace,[4 5])]);
padsize = matrix_size_complete(fps_to_xyz) - matrix_size_acquired(fps_to_xyz);
nTimes = size(kspace,4);
nCoils = size(kspace,5);
for iTime = 1:nTimes
    for iCoil = 1:nCoils
        % Zeropad, then ifftshift kspace
        padded_kspace = padarray(kspace(:,:,:,iTime,iCoil), 0.5*padsize, 0);
        padded_kspace = ifftshift(padded_kspace);

        % Take IFFT, then fftshift
        twistImage(:,:,:,iTime,iCoil) = fftshift(ifftn(padded_kspace));
    end
end

%% --- 8. Resolving Oversampling
% crop_amount = matrix_size_acquired-IMmatrix_crop_size;
% margin = floor(crop_amount ./ 2); 
% 
% final_IMspace = unTWISTed_IMspace(...
%     margin(1)+1 : margin(1)+IMmatrix_crop_size(1), ...
%     margin(2)+1 : margin(2)+IMmatrix_crop_size(2), ...
%     margin(3)+1 : margin(3)+IMmatrix_crop_size(3),...
%     :);


%% display
imslice(abs(twistImage))

% figure()
% subplot(1,3,1)
% h1 = imagesc(rot90(squeeze(abs(img_viaKspace(:,:,round(0.5*size(img_viaKspace,3)))))));
% colormap(gray)
% xlabel('R/L');
% ylabel('A/P');
% ax1 = ancestor(h1,'axes');
% set(ax1,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
%     'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([1 2 3]));
% subplot(1,3,2)
% h2 = imagesc(rot90(squeeze(abs(img_viaKspace(:,round(0.5*size(img_viaKspace,2)),:)))));
% colormap(gray)
% xlabel('R/L');
% ylabel('S/I');
% ax2 = ancestor(h2,'axes');
% set(ax2,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
%     'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([1 3 2]));
% title(encodingFullStr)
% subplot(1,3,3)
% h3 = imagesc(fliplr(rot90(squeeze(abs(img_viaKspace(round(0.5*size(img_viaKspace,1)),:,:))))));
% colormap(gray)
% xlabel('A/P');
% ylabel('S/I');
% ax3 = ancestor(h3,'axes');
% set(ax3,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
%     'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([2 3 1]));
% set(gcf,'Position',[680         665        1137         213]);
