clear all; 
close all; 
clc; 

%% 1) FOV and matrix size (scanner-style inputs)
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I 
encodingFullStr = formatEncodingString(freq_phase_slice);
disp(encodingFullStr)



load('Breast_Ultrafast_scan_parameters.mat') 
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

kOrderedIdx = [Sampling_Table.Frequency, Sampling_Table.("Row (phase)"), Sampling_Table.("Column (slice)")];

%% 3) Build WORLD k-space grid and map to the TWIST ordering
disp('Building WORLD k-space')
[kfreq_vec, kphase_vec, kslice_vec, ~, ~, ~] = ...
    computeKspaceGrid3D(FOV_acquired, matrix_size_acquired);
kfreq_orderedIdx = kOrderedIdx(:, 1);
kphase_orderedIdx = kOrderedIdx(:, 2);
kslice_orderedIdx = kOrderedIdx(:, 3);

ordKspace_freq = kfreq_vec(kfreq_orderedIdx);
ordKspace_phase = kphase_vec(kphase_orderedIdx);
ordKspace_slice = kslice_vec(kslice_orderedIdx);

k_fps = [ordKspace_freq, ordKspace_phase, ordKspace_slice]';
[k_xyz, fps_to_xyz] = mapKspaceFpsToXyz(k_fps, freq_phase_slice);

%% 5) Construct the breast phantom with the embedded enhancing vessel
t_PE = TR-(matrix_size_acquired(1)*dt_s); %TODO: rename these variables for clarity
t_s = t_PE+dt_s:dt_s:TR;
multipliers = 1:height(Sampling_Table)/matrix_size_acquired(1);

matrix_result = t_s(:) + (multipliers * TR);

Sampling_Table.Timing = matrix_result(:);

% % Force Timing to be increasing
% [sortedT,sortIdx] = sort(Sampling_Table.Timing);
% phantom = BreastPhantom(sortedT);

phantom = BreastPhantom(Sampling_Table.Timing);

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');

K_ordered = phantom.kspace(k_xyz(1, :), k_xyz(2, :), k_xyz(3, :));

Sampling_Table.Kspace_Value = K_ordered';

TWISTed_Kspace = zeros([matrix_size_acquired,max(Sampling_Table.Bj)+1]);
TWISTed_Kspace(Sampling_Table.("Linear Index")) = Sampling_Table.Kspace_Value;
%% ---  6. Updating K-Space from each measurement to undo TWIST

fprintf('\nUndoing TWIST...')

% Initialize the output array
unTWISTed_Kspace = zeros(size(TWISTed_Kspace));

% The first timepoint is the baseline for all coils
unTWISTed_Kspace(:,:,:,1,:) = TWISTed_Kspace(:,:,:,1,:);

% Loop only through timepoints (vectorized over coils)
for ii_timepoint = 2:size(TWISTed_Kspace,4)
    % 1. Create a temporary variable for the current timepoint's slice,
    % starting with data from the previous timepoint.
    dest_slice = unTWISTed_Kspace(:,:,:,ii_timepoint-1,:);

    % 2. Get the new sparse measurements for the current timepoint
    source_slice = TWISTed_Kspace(:,:,:,ii_timepoint,:);

    % 3. Create a logical mask of where the new measurements exist
    update_mask = (source_slice ~= 0);

    % 4. Use the mask to update the destination slice with the new values.
    dest_slice(update_mask) = source_slice(update_mask);

    % 5. Assign the updated slice back into the main array.
    unTWISTed_Kspace(:,:,:,ii_timepoint,:) = dest_slice;
end

clear('dest_slice','source_slice','update_mask')

%% --- 7. Resolving %resolution
padsize = matrix_size_complete - matrix_size_acquired;
final_kspace = padarray(unTWISTed_Kspace,padsize,0);

% Perform the inverse 3D FFT on all timepoints and coils at once
final_kspace_xyz = permute(final_kspace, [fps_to_xyz, 4]);
final_kspace_shifted = ifftshift(ifftshift(ifftshift(final_kspace_xyz, 1), 2), 3);
unTWISTed_IMspace = ifft(ifft(ifft(final_kspace_shifted, [], 1), [], 2), [], 3);
unTWISTed_IMspace = fftshift(fftshift(fftshift(unTWISTed_IMspace, 1), 2), 3);

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
imslice(abs(unTWISTed_IMspace))

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
