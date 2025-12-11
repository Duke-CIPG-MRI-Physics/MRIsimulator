%% DEMO: Analytic 3D MRI Simulation of basic shapes
%  ------------------------------------------------------
%  This demo shows how to:
%    1) Build a 3D k-space grid (WORLD coords, cycles/mm)
%    2) Evaluate the analytic Fourier transform of a various basic shapes
%    3) Reconstruct the 3D image using inverse FFT
%    4) Render the shape directly in IMAGE space using
%          estimateImage(xMesh,yMesh,zMesh)
%    5) Compare axial + coronal slices:
%          (A) FFT–reconstructed
%          (B) Geometry-based “exact shape” rendering
%
function demo_analtyicalShapes3D_temporalShapes()

%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [300 300 300];               % 300 mm cube
N      = [150 150 80];                % [Nx Ny Nz]

resolution = FOV_mm./N;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Build WORLD k-space grid
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

%% 3) Build centered WORLD image-domain meshgrids
[x_vec, y_vec, z_vec, ~, ~, ~] = computeCenteredImageGrid3D(FOV_mm, N);
midVol = round(N/2);

%% 4) Compute image coordinates for estimateImage()
fprintf('Precalculating pixels for estimateImage()...\n');
[x_ax, y_ax] = ndgrid(x_vec, y_vec);
z_ax  = z_vec(midVol(3)) * ones(size(x_ax));

[x_cor, z_cor] = ndgrid(x_vec, z_vec);
y_cor = y_vec(midVol(2)) * ones(size(x_cor));

%% 5) Instantiate analytic box (WORLD coordinates)
center = [0 0 0];
roll_deg   = 0;
pitch_deg  = 15;
yaw_deg    = 25;

%% 6) Construct array of shapes
center = [0 0 0];
roll_pitch_yaw = [0, 15, 35];
intensity = 1;

t_samples = reshape(1:numel(kx),size(kx));

% boxParams = struct('Lx_mm', 60, 'Ly_mm', 40, 'Lz_mm', 20); % [mm];
cylParams = struct('radius_mm', 60, 'length_mm', 80);
ellipParams = struct('a_mm', 140, 'b_mm', 100, 'c_mm', 70);
ellipCylParams = struct('a_mm', 140, 'b_mm', 100, 'length_mm', 70);
sphereParams = struct('radius_mm', 60);

shapes = [AnalyticalBox3D(@()calcBoxParams(t_samples), intensity, center, roll_pitch_yaw),...
    AnalyticalCylinder3D(cylParams, intensity, center, roll_pitch_yaw),...
    AnalyticalEllipsoid3D(ellipParams, intensity, center, roll_pitch_yaw),...
    AnalyticalEllipticalCylinder3D(ellipCylParams, intensity, center, roll_pitch_yaw),...
    AnalyticalSphere3D(sphereParams, intensity, center, roll_pitch_yaw)];

%% 7) Loop through shapes
nShapes = length(shapes);
for iShape = 1:nShapes
    % 7a) Compute class name
    thisShape = shapes(iShape);
    thisClassName = class(thisShape);

    % 7a) Compute analytic k-space
    fprintf('Evaluating analytic k-space...\n');
    K = thisShape.kspace(kx, ky, kz);


    % 7b) Reconstruct 3D image via inverse FFT
    fprintf('Performing inverse FFT...\n');
    img_viaKspace = fftshift(ifftn(ifftshift(K)));

    % 8) Visualization: *four* images
    figure('Name',[thisClassName ': k-space rendered over time'],'Color','w');

    % --- Axial, FFT ---
    subplot(1,2,1);
    imagesc(x_vec, y_vec, abs(squeeze(img_viaKspace(:,:,midVol(3)))).');
    axis image; colormap gray; colorbar;
    set(gca,'YDir','normal');
    xlabel('x (mm)'); ylabel('y (mm)');
    title('Axial (k-space rendered)');

    % --- Coronal, FFT ---
    subplot(1,2,2);
    imagesc(y_vec, z_vec, abs(squeeze(img_viaKspace(midVol(1),:,:))).');
    axis image; colormap gray; colorbar;
    set(gca,'YDir','normal');
    xlabel('y (mm)'); ylabel('z (mm)');
    title('Coronal (k-space rendered)');

    sgtitle([thisClassName ': k-space rendered over time'],'FontSize',16);

end



    function boxParams = calcBoxParams(t)
        boxParams = struct('Lx_mm', 60 - 10*t/max(t(:)), 'Ly_mm', 40+70*t/max(t(:)), 'Lz_mm', 20*ones(size(t))); % [mm];
    end
end