% 3D phantom: single cylinder using analytic k-space.

clear; clc;

%% FOV and sampling
FOV_mm = 300;     % 300 mm

dx_mm = 2; 
dy_mm = 2; 
dz_mm = 6;   

Nx = FOV_mm / dx_mm;      % 300
Ny = FOV_mm / dy_mm;      % 300
Nz = FOV_mm / dz_mm;      % 300

if any(mod([Nx Ny Nz],1) ~= 0)
    error('FOV must be divisible by voxel size to get integer grid sizes.');
end

fprintf('Grid size: %d x %d x %d\n', Nx, Ny, Nz);

%% Spatial frequency grid (cycles/mm)
kx_vec = computeKspaceAxis(FOV_mm, Nx);
ky_vec = computeKspaceAxis(FOV_mm, Ny);
kz_vec = computeKspaceAxis(FOV_mm, Nz);

[kx, ky, kz] = ndgrid(kx_vec, ky_vec, kz_vec);

%% Define cylinder geometry (body axis along z)
R_mm = 20;      % radius
L_mm = 100;     % length along z

center_cyl  = [0 0 0];   % world center (mm)
roll_deg  = 0;
pitch_deg = 0;
yaw_deg   = 0;

% Requires AnalyticalCylinder3D + AnalyticalShape3D on path
cyl = AnalyticalCylinder3D(center_cyl, ...
                                        roll_deg, pitch_deg, yaw_deg, ...
                                        R_mm, L_mm);

%% Compute analytic k-space for cylinder
fprintf('Evaluating analytic k-space (cylinder)...\n');
K_total = cyl.kspace(kx, ky, kz);

%% 3D inverse FFT to obtain centered image
fprintf('Performing 3D iFFT...\n');

img_complex = fftshift(ifftn(ifftshift(K_total)));
img_mag     = abs(img_complex);
img_mag     = img_mag / max(img_mag(:));

%% Spatial coordinate grid (mm), centered around 0
x_vec = ((0:Nx-1) - floor(Nx/2)) * dx_mm;  % ~[-FOV/2, +FOV/2)
y_vec = ((0:Ny-1) - floor(Ny/2)) * dy_mm;
z_vec = ((0:Nz-1) - floor(Nz/2)) * dz_mm;

%% Display central slices
midZ = round(Nz/2);   % z ≈ 0
midY = round(Ny/2);   % y ≈ 0
midX = round(Nx/2);   % x ≈ 0



global_clim = [0 1];   % or [0 max(img_mag(:))], but you've normalized

figure;
for iS = 1:Nz
    subplot(ceil(sqrt(Nz)), ceil(sqrt(Nz)), iS);
    imagesc(squeeze(img_mag(:,:,iS)).');
    colormap gray;
    caxis(global_clim);        % <-- critical: same scale for all slices
    axis image off;
end

% figure;
% for iS = 1:Nz
% subplot(ceil(sqrt(Nz)),ceil(sqrt(Nz)),iS);
% imagesc( squeeze(img_mag(:,:,iS)).');
% colormap gray;
% end

figure;
subplot(1,2,1);
imagesc(x_vec, z_vec, squeeze(img_mag(:,midY,:)).');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('x (mm)'); ylabel('z (mm)');
title('Sagittal (y \approx 0)');

subplot(1,2,2);
imagesc(y_vec, z_vec, squeeze(img_mag(midX,:,:)).');
axis image; colormap gray; colorbar;
set(gca,'YDir','normal');
xlabel('y (mm)'); ylabel('z (mm)');
title('Coronal (x \approx 0)');