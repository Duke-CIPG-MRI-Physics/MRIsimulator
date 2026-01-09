% Start with a clean workspace
clc; clear all; close all

% Load example volume
load("xcat_MRI_3D.mat");

% Simulate coils
nCoilsEven = 4;
sigma = [0.8 0.55];
coilSplitDirection = 'z';
coilDirection = 'y';
opts = struct;
opts.planeOffset = 1.02;
opts.phaseAmp = 0.6;
opts.noiseRMS = 1e-2;     % adjust for desired SNR
opts.dtype           = "single";      % or "double"    (must be *string*)

[Combined_output_coils, Combined_output_coils_kspace, coilSens] = ...
    simulateCoils(Combined_output, nCoilsEven, sigma, ...
    coilSplitDirection, coilDirection);

figure
imshow3D(squeeze(abs(Combined_output_coils(:,:,:,1))))
figure
imshow3D(angle(Combined_output_coils(:,:,:,1)))