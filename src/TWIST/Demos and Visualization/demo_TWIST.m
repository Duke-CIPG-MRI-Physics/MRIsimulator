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


[A1_idx_outerIn, A2_idx_innerOut] = calculateTwistSamplingOrder2D(pA, pB, FOV_acquired, matrix_size_acquired, R);
