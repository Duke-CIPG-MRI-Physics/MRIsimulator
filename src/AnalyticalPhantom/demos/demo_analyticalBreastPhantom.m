%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [500 500 400];
N      = [150 150 150]; % [Nx Ny Nz]
resolution = FOV_mm./N;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Configure time-varying contrast wash-in for the embedded vessel
t_s = linspace(0, 5, 25).';
vesselRadius_mm = 2.5;
total_vessel_length_mm = 100;

% Linear contrast volume so the vessel is fully enhanced at the final time
totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
V_contrast_mm3 = linspace(0, totalVolume_mm3, numel(t_s)).';

enhancedLength_mm = computeContrastWashIn(t_s, vesselRadius_mm, V_contrast_mm3);

%% 3) Construct the breast phantom with the embedded enhancing vessel
phantom = BreastPhantom(t_s, V_contrast_mm3, vesselRadius_mm);

%% 4) Build WORLD k-space grid
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

%% 5) Compute analytic k-space for the phantom
fprintf('Evaluating analytic k-space...\n');
K = phantom.kspace(kx, ky, kz);

%% 6) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% display
imslice(abs(img_viaKspace))

figure;
plot(t_s, min(enhancedLength_mm, total_vessel_length_mm), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Enhanced vessel length [mm]');
title('Linear contrast wash-in to full vessel length');
grid on;
