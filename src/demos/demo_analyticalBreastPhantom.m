clear all; 
close all; 
clc; 
restoredefaultpath; 
run('C:\code\mriSimulator\cipg_setup.m');
run('C:\code\Duke-CIVM-MRI-Tools\setup.m');
savepath

%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [400 400 400];
N      = [200 200 30]; % [Nx Ny Nz]
resolution = FOV_mm./N;
resolution_normalized = resolution/max(resolution);

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Configure acquisition ordering and timing
freq_phase_slice = [3 2 1]; % 1 = R/L, 2=A/P, 3 = S/I
[kOrderedIdx, ~] = orderRectilinearKspace(N, freq_phase_slice, 1, 1);

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

%% 3) Build WORLD k-space grid and map to the rectilinear ordering
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);
kx_orderedIdx = kOrderedIdx(1, :);
ky_orderedIdx = kOrderedIdx(2, :);
kz_orderedIdx = kOrderedIdx(3, :);

ordKsx_kx = kx_vec(kx_orderedIdx)';
ordKsx_ky = ky_vec(ky_orderedIdx)'; 
ordKsx_kz = kz_vec(kz_orderedIdx)';

%% 5) Construct the breast phantom with the embedded enhancing vessel
phantom = BreastPhantom(0);

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');
K_ordered = phantom.kspace(ordKsx_kx, ordKsx_ky, ordKsx_kz);

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
h2 = imagesc(rot90(squeeze(abs(img_viaKspace(:,116,:)))));
colormap(gray)
xlabel('R/L');
ylabel('S/I');
ax2 = ancestor(h2,'axes');
set(ax2,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
    'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([1 3 2]));
title(encodingFullStr)
subplot(1,3,3)
h3 = imagesc(fliplr(rot90(squeeze(abs(img_viaKspace(106,:,:))))));
colormap(gray)
xlabel('A/P');
ylabel('S/I');
ax3 = ancestor(h3,'axes');
set(ax3,'XTick',[],'XTickLabel',[], 'YTick',[],'YTickLabel',[],...
    'PlotBoxAspectRatioMode','auto','DataAspectRatio',resolution_normalized([2 3 1]));
set(gcf,'Position',[680         665        1137         213]);
