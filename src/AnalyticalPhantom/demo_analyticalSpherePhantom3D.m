%% DEMO: Analytic 3D MRI Simulation of a Solid Sphere
%  --------------------------------------------------
%  This demo shows how to:
%    1) Build a 3D k-space grid (WORLD coords, cycles/mm)
%    2) Evaluate the analytic Fourier transform of a sphere
%    3) Reconstruct the 3D image using inverse FFT
%    4) Render the sphere directly in IMAGE space using
%          estimateImageShape(xMesh,yMesh,zMesh)
%    5) Compare axial + coronal slices:
%          (A) FFT–reconstructed
%          (B) Geometry-based “exact shape” rendering
%
clear; clc;

%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [300 300 300];               % 300 mm cube
N      = [150 150 80];      % [Nx Ny Nz]

resolution = FOV_mm./N;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Build WORLD k-space grid
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

%% 3) Instantiate analytic sphere (WORLD coordinates)
R_mm = 140;
center_sph = [0 0 0];

sph = AnalyticalSphere3D(R_mm, [], center_sph, [0 0 0]);

%% 4) Compute analytic k-space
fprintf('Evaluating analytic k-space of sphere...\n');
K = sph.kspace(kx, ky, kz);

%% 5) Reconstruct 3D image via inverse FFT
fprintf('Performing inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% 6) Build centered WORLD image-domain meshgrids
[x_vec, y_vec, z_vec, xMesh, yMesh, zMesh] = computeCenteredImageGrid3D(FOV_mm, N);
midVol = round(N/2);

%% 7) Compute geometry-based “exact” image using percentInsideShape()
fprintf('Rasterizing shape using estimateImageShape()...\n');

% Extract 2D axial and coronal slices in WORLD coordinates
[x_ax, y_ax] = ndgrid(x_vec, y_vec);
z_ax  = z_vec(midVol(3)) * ones(size(x_ax));

[x_cor, z_cor] = ndgrid(x_vec, z_vec);
y_cor = y_vec(midVol(2)) * ones(size(x_cor));

% Get fractional volumes (0–1)
frac_ax  = sph.estimateImageShape(x_ax, y_ax, z_ax);
frac_cor = sph.estimateImageShape(x_cor, y_cor, z_cor);

% Normalize (optional, matches FFT-range)
frac_ax  = frac_ax / max(frac_ax(:));
frac_cor = frac_cor / max(frac_cor(:));

%% 8) Visualization: *four* images
figure('Name','Sphere: FFT vs Image-Space Shape Rendering','Color','w');

% --- Axial, FFT ---
subplot(2,2,1);
imagesc(x_vec, y_vec, abs(squeeze(img_viaKspace(:,:,midVol(3)))).');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('x (mm)'); ylabel('y (mm)');
title('Axial (FFT Reconstructed)');

% --- Axial, Shape-rendered ---
subplot(2,2,2);
imagesc(x_vec, y_vec, frac_ax.');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('x (mm)'); ylabel('y (mm)');
title('Axial (Direct Shape Rendering)');

% --- Coronal, FFT ---
subplot(2,2,3);
imagesc(y_vec, z_vec, abs(squeeze(img_viaKspace(midVol(1),:,:))).');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('y (mm)'); ylabel('z (mm)');
title('Coronal (FFT Reconstructed)');

% --- Coronal, Shape-rendered ---
subplot(2,2,4);
imagesc(y_vec, z_vec, frac_cor.');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('y (mm)'); ylabel('z (mm)');
title('Coronal (Direct Shape Rendering)');

sgtitle('Analytic Sphere: k-space FFT vs Geometry-based Rendering','FontSize',16);
