clear all; 
close all; 
clc; 
restoredefaultpath; 
run('C:\code\mriSimulator\cipg_setup.m');
run('C:\code\Duke-CIVM-MRI-Tools\setup.m');
savepath

%% 1) FOV and matrix size (scanner-style inputs)
FOV_mm = [500 500 400];
N      = [150 150 150]; % [Nx Ny Nz]
resolution = FOV_mm./N;

fprintf('Voxel size: [%.2f %.2f %.2f] mm\n', resolution(1), resolution(2), resolution(3));
fprintf('Grid size : %d x %d x %d\n', N(1), N(2), N(3));

%% 2) Configure acquisition ordering and timing
freq_phase_slice = [2 1 3];
dt = 4e-6;   % dwell time between frequency-encode samples [s]
TR = 5e-3;   % time between starts of successive frequency-encode lines [s]

[kOrderedIdx, tSamp] = orderRectilinearKspace(N, freq_phase_slice, dt, TR);
t_s = tSamp(:); % use sampling timestamps as the phantom time base

%% 3) Build WORLD k-space grid and map to the rectilinear ordering
[kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N);
kx_orderedIdx = kOrderedIdx(1, :);
ky_orderedIdx = kOrderedIdx(2, :);
kz_orderedIdx = kOrderedIdx(3, :);

%% 4) Configure time-varying contrast wash-in for the embedded vessel
vesselRadius_mm = 2.5;
total_vessel_length_mm = 100;
totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;

startTime = 0.25 * t_s(end);
endTime   = 0.75 * t_s(end);

V_contrast_mm3 = zeros(numel(t_s), 1);
midRamp = t_s >= startTime & t_s <= endTime;
V_contrast_mm3(midRamp) = totalVolume_mm3 * (t_s(midRamp) - startTime) ./ (endTime - startTime);
V_contrast_mm3(t_s > endTime) = totalVolume_mm3;

enhancedLength_mm = computeContrastWashIn(t_s, vesselRadius_mm, V_contrast_mm3);

%% 5) Construct the breast phantom with the embedded enhancing vessel
phantom = BreastPhantom(t_s, V_contrast_mm3, vesselRadius_mm);

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');
K_ordered = phantom.kspace(kx_vec(kx_orderedIdx)', ky_vec(ky_orderedIdx)', kz_vec(kz_orderedIdx)');

% Reassemble onto the kx/ky/kz grid for reconstruction
K = zeros(N);
linearIdx = sub2ind(N, kOrderedIdx(1, :), kOrderedIdx(2, :), kOrderedIdx(3, :));
K(linearIdx) = K_ordered;

%% 7) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% display
imslice(abs(img_viaKspace))

figure;
plot(t_s, min(enhancedLength_mm, total_vessel_length_mm), '-','LineWidth', 2);
hold on
plot(t_s, min(enhancedLength_mm, total_vessel_length_mm), '.r');
xlabel('Time [s]');
ylabel('Enhanced vessel length [mm]');
title('Linear contrast wash-in to full vessel length');
grid on;
