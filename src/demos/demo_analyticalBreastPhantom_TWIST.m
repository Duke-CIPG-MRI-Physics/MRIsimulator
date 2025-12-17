clear all; 
close all; 
clc; 
% restoredefaultpath; 
% run('C:\code\mriSimulator\cipg_setup.m');
% run('C:\code\Duke-CIVM-MRI-Tools\setup.m');
% savepath

%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [400 400 400];
N      = [100 100 100]; % [Nx Ny Nz]
resolution = FOV_mm./N;
resolution_normalized = resolution/max(resolution);

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Configure acquisition ordering and timing
freq_phase_slice = [3 2 1]; % 1 = R/L, 2=A/P, 3 = S/I
% [kOrderedIdx, ~] = orderRectilinearKspace(N, freq_phase_slice, 1, 1);

% parameters from clinical ultrafast protocol
FOV_read_mm = 350;
FOV_phase_pct = 100;
oversampling_phase_pct = 20; 
oversampling_slice_pct = 33.3; 
slices_per_slab = 240;
slice_thickness_mm = 1; % Note this is the reconstructed slice thickness, not the nominal slice thickness 
base_resolution = 224;
phase_resolution_pct = 100;
slice_resolution_pct = 80;
freq_phase_slice = [3 2 1]; % 1 = R/L, 2=A/P, 3 = S/I

% TODO - interpolation is off in current protocol, but we could consider it as an option in the future...

% Contrast parameters
rBW_HzPerPix = 570;
TR = 5.88E-3;   
TE = 2.63E-3;

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*base_resolution;
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]


pA = 0.05;
Nb = 10;
Time_Measured = 90; %sec
R = 1;
PF_Factor = 1;

Sampling_Table = Ultrafast_Sampling(N,pA,Nb,Time_Measured,TR,R,PF_Factor);

kOrderedIdx = [Sampling_Table.Frequency, Sampling_Table.("Row (phase)"), Sampling_Table.("Column (slice)")];


encodingDimStr = {'freq:',', phase:',', slice:'};
encodingFullStr = '';
for iDim = 1:3
    switch freq_phase_slice(iDim)
        case 1
            thisDimStr = 'R/L';
        case 2
            thisDimStr = 'A/P';
        case 3
            thisDimStr = 'S/I';
        otherwise
            error('freq_phase_slice can only be numbers 1-3')
    end
    encodingFullStr = [encodingFullStr encodingDimStr{iDim} thisDimStr];
end
encodingFullStr

%% 3) Build WORLD k-space grid and map to the TWIST ordering
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

kx_orderedIdx = kOrderedIdx(:, 1);
ky_orderedIdx = kOrderedIdx(:, 2);
kz_orderedIdx = kOrderedIdx(:, 3);

ordKsx_kx = kx_vec(kx_orderedIdx);
ordKsx_ky = ky_vec(ky_orderedIdx); 
ordKsx_kz = kz_vec(kz_orderedIdx);

%% 5) Construct the breast phantom with the embedded enhancing vessel
t_PE = TR-(N(1)*dt_s);
t_s = t_PE+dt_s:dt_s:TR;
multipliers = 1:height(Sampling_Table)/N(1);

matrix_result = t_s(:) + (multipliers * TR);

Sampling_Table.Timing = matrix_result(:);

phantom = BreastPhantom(Sampling_Table.Timing);

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');
K_ordered = phantom.kspace(ordKsx_kx, ordKsx_ky, ordKsx_kz);

% TODO1 - fix the recon
% Either break this into separate time-recon windows 

% Reassemble onto the kx/ky/kz grid for reconstruction
K = zeros(N);
linearIdx = sub2ind(N, kx_orderedIdx, ky_orderedIdx, kz_orderedIdx);
K(linearIdx) = K_ordered;

%% 6) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));



%% display
imslice(abs(img_viaKspace))

figure()
subplot(1,3,1)
h1 = imagesc(rot90(squeeze(abs(img_viaKspace(:,:,round(0.5*size(img_viaKspace,3)))))));
colormap(gray)
xlabel('R/L');
ylabel('A/P');
ax1 = ancestor(h1,'axes');
set(ax1,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
    'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([1 2 3]));
subplot(1,3,2)
h2 = imagesc(rot90(squeeze(abs(img_viaKspace(:,round(0.5*size(img_viaKspace,2)),:)))));
colormap(gray)
xlabel('R/L');
ylabel('S/I');
ax2 = ancestor(h2,'axes');
set(ax2,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
    'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([1 3 2]));
title(encodingFullStr)
subplot(1,3,3)
h3 = imagesc(fliplr(rot90(squeeze(abs(img_viaKspace(round(0.5*size(img_viaKspace,1)),:,:))))));
colormap(gray)
xlabel('A/P');
ylabel('S/I');
ax3 = ancestor(h3,'axes');
set(ax3,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
    'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([2 3 1]));
set(gcf,'Position',[680         665        1137         213]);
