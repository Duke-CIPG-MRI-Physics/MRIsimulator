%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [500 500 400];
N      = [150 150 150]; % [Nx Ny Nz]
resolution = FOV_mm./N;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% Define starting parameters
bodyShift = -80;

% Heart
heart_center = [0 bodyShift 0]; % mm, anywhere you want in WORLD coords
heart_a_mm = 50;
heart_b_mm = 27;
heart_c_mm = 30;

heart_roll_deg  = 0;% does nothing
heart_pitch_deg = -65;
heart_yaw_deg   = 70;

heart = AnalyticalEllipsoid3D(heart_center, ...
    heart_roll_deg, heart_pitch_deg, heart_yaw_deg, ...
    heart_a_mm, heart_b_mm, heart_c_mm);

% Left lung
lung_a_mm = 140;
lung_b_mm = 50;
lung_c_mm = 50;

lung_roll_R  = 0;
lung_pitch_R = 95;
lung_yaw_R   = 0; % does nothing

lungSeparation = max(lung_b_mm,lung_c_mm) + max(heart_b_mm,heart_c_mm) + 2;
center_R_mm = [lungSeparation, bodyShift, 0];

rightLung = AnalyticalEllipsoid3D(center_R_mm, ...
    lung_roll_R, lung_pitch_R, lung_yaw_R, ...
    lung_a_mm, lung_b_mm, lung_c_mm);

lung_roll_L  = 0;
lung_pitch_L = 85;
lung_yaw_L   = 0;

center_L_mm = [-lungSeparation, bodyShift, 0];

leftLung = AnalyticalEllipsoid3D(center_L_mm, ...
    lung_roll_L, lung_pitch_L, lung_yaw_L, ...
    lung_a_mm, lung_b_mm, lung_c_mm);

% Peripheral Fat
bodyCenter = [0 bodyShift 0];
fat_roll_deg = 0;
fat_pitch_deg = 0;
fat_yaw_deg = 0;
fatThickness_mm = 10;
tissueThickness_mm = 10;

patientThickness_outer_mm = 1.85*(max(lung_b_mm, lung_c_mm)+fatThickness_mm);
patientWidth_outer_mm = 1.5*(lungSeparation+max(lung_b_mm, lung_c_mm)+fatThickness_mm+tissueThickness_mm )+2;

fat_outer = AnalyticalEllipticalCylinder3D(bodyCenter, fat_roll_deg, fat_pitch_deg, fat_yaw_deg, ...
    patientWidth_outer_mm, patientThickness_outer_mm, 0.9*FOV_mm(3));

patientThickness_inner_mm = patientThickness_outer_mm - 2*fatThickness_mm;
patientWidth_inner_mm = patientWidth_outer_mm - 2*fatThickness_mm;

fat_inner = AnalyticalEllipticalCylinder3D(bodyCenter, fat_roll_deg, fat_pitch_deg, fat_yaw_deg, ...
    patientWidth_inner_mm, patientThickness_inner_mm, 0.9*FOV_mm(3));

% right breast
breast_gap_mm = 50;
breast_roll_deg = 0;
breast_pitch_deg = 90;
breast_yaw_deg = 90;

breast_radius_mm = 65;
breast_depth_mm = 200;

right_breast_center = [breast_radius_mm+0.5*breast_gap_mm, ...
    bodyShift + 0.5*(breast_depth_mm) + patientThickness_outer_mm, ...
    0];

breast_right = AnalyticalCylinder3D(right_breast_center, breast_roll_deg, breast_pitch_deg, breast_yaw_deg, ...
    breast_radius_mm, breast_depth_mm);

% left breast
left_breast_center = [-right_breast_center(1) right_breast_center(2:3)];
breast_left = AnalyticalCylinder3D(left_breast_center, breast_roll_deg, breast_pitch_deg, breast_yaw_deg, ...
    breast_radius_mm, breast_depth_mm);


%% 2) Build WORLD k-space grid
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);

%% 4) Compute analytic k-space for the cylinder
fprintf('Evaluating analytic k-space...\n');
fatIntensity = 2;
tissueIntensity = 0.5;
heartIntensity = 1;
lungIntensity = 0.1;
breast_intensity = 0.5;
K = fatIntensity*(fat_outer.kspace(kx, ky, kz) - fat_inner.kspace(kx, ky, kz)) + ...
    tissueIntensity*(fat_inner.kspace(kx, ky, kz) - ...
    (heart.kspace(kx, ky, kz) + rightLung.kspace(kx, ky, kz) + leftLung.kspace(kx, ky, kz))) + ...
    heartIntensity*heart.kspace(kx, ky, kz) + ...
    lungIntensity*rightLung.kspace(kx, ky, kz) + ...
    lungIntensity*leftLung.kspace(kx, ky, kz) + ...
    breast_intensity*breast_left.kspace(kx, ky, kz) + ...
    breast_intensity*breast_right.kspace(kx, ky, kz);

%% 5) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% display
imslice(abs(img_viaKspace))