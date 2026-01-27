clc; clear all; close all; 

freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I 
encodingFullStr = formatEncodingString(freq_phase_slice);
disp(encodingFullStr)

load('Breast_Ultrafast_scan_parameters.mat')
[FOV_acquired,matrix_size_complete,matrix_size_acquired,voxel_size_mm,nyquist_resolution_mm,IMmatrix_crop_size] = convert_Siemens_parameters(scan_parameters);




% Contrast parameters
rBW_HzPerPix = 570;
TR = (5.88E-3);  
TE = 2.63E-3;

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*matrix_size_acquired(1);  
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

pA = 0.05;
pB = 0.1;

R = [2 3];
PF_Factor = [6/8 6/8];


[A1_idx_outerIn, A2_idx_innerOut, B1_idx_innerOut,B2_idx_outerIn] = ...
    calculateTwistSamplingOrder2D(pA, pB, FOV_acquired, matrix_size_acquired, R, PF_Factor);

%% --- Phase/Slice linear indices for downstream use
phase_slice_size = matrix_size_acquired(2:3);
twist_phase_slice_indices = struct( ...
    'Size', phase_slice_size, ...
    'A1_outerIn', A1_idx_outerIn, ...
    'A2_innerOut', A2_idx_innerOut, ...
    'B1_innerOut', {B1_idx_innerOut}, ...
    'B2_outerIn', {B2_idx_outerIn});

%% --- Frequency-encoded linear indices
n_freq = matrix_size_acquired(1);
freq_offsets = (1:n_freq)';
full_phase_slice_size = [n_freq phase_slice_size];
linearize_freq = @(idx) reshape((idx(:) - 1) * n_freq + freq_offsets, [], 1);

twist_frequency_indices = struct( ...
    'Size', full_phase_slice_size, ...
    'A1_outerIn', linearize_freq(A1_idx_outerIn), ...
    'A2_innerOut', linearize_freq(A2_idx_innerOut), ...
    'B1_innerOut', {cellfun(linearize_freq, B1_idx_innerOut, 'UniformOutput', false)}, ...
    'B2_outerIn', {cellfun(linearize_freq, B2_idx_outerIn, 'UniformOutput', false)});
