%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [500 500 400];
N      = [150 150 150]; % [Nx Ny Nz]
resolution = FOV_mm./N;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Construct the lung phantom
phantom = LungPhantom();

%% 3) Build WORLD k-space grid
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

%% 4) Compute analytic k-space for the phantom
fprintf('Evaluating analytic k-space...\n');
K = phantom.kspace(kx, ky, kz);

%% 5) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% display
imslice(abs(img_viaKspace))
