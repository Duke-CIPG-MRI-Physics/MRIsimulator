clear all; 
close all; 
clc; 

%% 1) FOV and matrix size (scanner-style inputs)
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I 
encodingFullStr = formatEncodingString(freq_phase_slice);
disp(encodingFullStr)

breastPhantomParams = createBreastPhantomParams();

load('fast_scan_parameters.mat') 
% load('Breast_Ultrafast_scan_parameters.mat')
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

% Calculate the time for a single twist frame
samplesInFrame = histcounts(Sampling_Table.Bj,-0.5:1:(max(Sampling_Table.Bj(:))+0.5));
samplesPerFrame = max(samplesInFrame(2:end)); % ignore first frame which measures all data
framesPerRecon = ceil(samplesInFrame(1)/samplesPerFrame);
trsPerFrame = samplesPerFrame/scan_parameters("Base Resolution");
timePerFrame = trsPerFrame*TR;
timePerRecon = framesPerRecon*timePerFrame;


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
clear TR_counts

% Update startInjectionTime_s to be relative to first frame ending
endOfSecondFrame = max(Sampling_Table(Sampling_Table.Bj == 0,:).Timing);

% Injected contrast parameters
breastPhantomParams.startEnhancement_s = breastPhantomParams.startInjectionTime_s + endOfSecondFrame;
breastPhantomParams.enhancementDuration_s = 10;
breastPhantomParams.unenhancedIntensity = 0.4;
breastPhantomParams.enhancedIntensity = 0.4;
breastPhantomParams.lesionIntensityFunction = @(t_s) min(2, max(0, 2 .* ...
    (t_s-breastPhantomParams.startEnhancement_s) ./ ...
    breastPhantomParams.enhancementDuration_s)) - ...
    breastPhantomParams.breastIntensity; % subtracting intensity allows us to avoid needing to subtract shapes

figure();
plot(Sampling_Table.Timing(:),breastPhantomParams.lesionIntensityFunction(Sampling_Table.Timing(:)) ...
    + breastPhantomParams.breastIntensity);


%% Contrast timing visualization

phantom = BreastPhantom(breastPhantomParams);

%% Perform TWIST, calculating two time frames at a time to minimize memory overhead
maxChumkSize = 500000;
previousMask = (Sampling_Table.Bj == 0);


nTimes = max(Sampling_Table.Bj)+1;
fprintf(['Reconstructing TWIST time %d of %d (%.1f%% complete).\n'], ...
        1, nTimes, 0/nTimes*100);
currentKspace = nan(matrix_size_acquired);
currentIdx = sub2ind(matrix_size_acquired, ...
    Sampling_Table.Frequency(previousMask), ...          
    Sampling_Table.("Row (phase)")(previousMask), ...    
    Sampling_Table.("Column (slice)")(previousMask));
currentKspace(currentIdx) = phantom.kspaceAtTime(k_spatFreq_xyz(1, previousMask), ...
    k_spatFreq_xyz(2, previousMask), ...
    k_spatFreq_xyz(3, previousMask), ...
    Sampling_Table.Timing(previousMask)', ...
    maxChumkSize)';

% initialize TWIST image with first frame by permuting to XYZ from FPS, 
% zeropading, ifftshifting k-space, taking IFFT, and fftshifting to get image
twistImage = zeros([matrix_size_complete(fps_to_xyz) nTimes]);
padsize = matrix_size_complete(fps_to_xyz) - matrix_size_acquired(fps_to_xyz);
twistImage(:,:,:,1) = fftshift(ifftn(ifftshift(padarray(...
        permute(currentKspace, [fps_to_xyz]), 0.5*padsize, 0)))); 
previousKspace = currentKspace;

for iTime = 2:nTimes
    fprintf(['Reconstructing TWIST time %d of %d (%.1f%% complete).\n'], ...
        iTime, nTimes, (iTime-1)/nTimes*100);

    % Calculate current Kspace Samples, putting k-space points in correct locations
    currentMask = (Sampling_Table.Bj == (iTime - 1));
    currentKspace = nan(matrix_size_acquired);

    currentIdx = sub2ind(matrix_size_acquired, ...
        Sampling_Table.Frequency(currentMask), ...
        Sampling_Table.("Row (phase)")(currentMask), ...
        Sampling_Table.("Column (slice)")(currentMask));
    currentKspace(currentIdx) = phantom.kspaceAtTime(k_spatFreq_xyz(1, currentMask), ...
        k_spatFreq_xyz(2, currentMask), ...
        k_spatFreq_xyz(3, currentMask), ...
        Sampling_Table.Timing(currentMask)', ...
        maxChumkSize)';
    
    % Fill in missing k-space points. 
    currentKspace(isnan(currentKspace)) = previousKspace(isnan(currentKspace));

    % Permute to XYZ from FPS, zeropad, then ifftshift k-space, then take
    % IFFT and fftshift to get image
    twistImage(:,:,:,iTime) = fftshift(ifftn(ifftshift(padarray(...
        permute(currentKspace, [fps_to_xyz]), 0.5*padsize, 0))));     

    % prepare for next TWIST frame
    previousKspace = currentKspace;

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
