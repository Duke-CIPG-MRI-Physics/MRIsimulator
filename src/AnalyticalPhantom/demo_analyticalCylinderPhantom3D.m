%% DEMO: Analytic 3D MRI Simulation of a Finite Cylinder
%  This demo illustrates:
%   - How to build a 3D k-space grid from FOV + matrix size
%   - How to use AnalyticalCylinder3D + AnalyticalShape3D
%   - How to reconstruct and visualize the resulting 3D phantom
%
% Things to notice
% - The cylinder is finite along z: only slices with |z| < L_mm/2 are bright.
% - Oscillatory sidelobes along z (and radially) are expected due to:
%       finite k-space sampling (i.e., band-limiting) → sinc/Gibbs ringing.
% - The object is never voxelized explicitly; we work entirely in k-space.
% - Multiple shapes can be added by summing their analytic FTs in K:
%       K_total = cyl1.kspace(...) + cyl2.kspace(...) + ...
clear; clc;

%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = 300;               % scalar → 300 mm in x,y,z
N      = [150 150 50];      % [Nx Ny Nz]; implies anisotropic voxels

Nx = N(1); Ny = N(2); Nz = N(3);

dx_mm = FOV_mm / Nx;
dy_mm = FOV_mm / Ny;
dz_mm = FOV_mm / Nz;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', dx_mm, dy_mm, dz_mm);
fprintf('Grid size : %d x %d x %d voxels\n', Nx, Ny, Nz);

%% 2) Build k-space grid (WORLD coordinates, cycles/mm)
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

%% 3) Define cylinder geometry (in WORLD coords)
% BODY frame: cylinder axis along +z, radius R_mm, length L_mm
% WORLD frame: specified by center + Euler angles

R_mm = 70;      % radius [mm]
L_mm = 100;     % length [mm]

center_cyl = [0 0 0];   % world center [mm]
roll_deg   = 0;         % rotations (deg)
pitch_deg  = 0;
yaw_deg    = 0;

% Create the analytic cylinder shape
cyl = AnalyticalCylinder3D(center_cyl, ...
                           roll_deg, pitch_deg, yaw_deg, ...
                           R_mm, L_mm);

% (Optional) Example of updating geometry:
% cyl.setRadius(30);
% cyl.setLength(120);
% cyl.setCenter([20 0 0]);
% cyl.setOrientation(0, 0, 45);

%% 4) Compute analytic k-space for the cylinder
fprintf('Evaluating analytic k-space...\n');
K = cyl.kspace(kx, ky, kz);

%% 5) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');

img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% 6) Build centered spatial grids (WORLD coords)
[x_vec, y_vec, z_vec, ~, ~, ~] = computeCenteredImageGrid3D(FOV_mm, N);

midX = round(Nx/2);
midY = round(Ny/2);
midZ = round(Nz/2);

%% 7) axial and coronal views
figure('Name','Sagittal & Coronal');

subplot(1,2,1);
imagesc(x_vec, y_vec, abs(squeeze(img_viaKspace(:,:,midZ))).');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('x (mm)');
ylabel('y (mm)');
title('axial slice (z \approx 0)');

subplot(1,2,2);
imagesc(y_vec, z_vec, abs(squeeze(img_viaKspace(midX,:,:))).');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('y (mm)');
ylabel('z (mm)');
title('Coronal slice (x \approx 0)');
