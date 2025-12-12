%% demo_breathing_motion.m
% End-to-end demo: define breathing params, compute geometry, visualize.

clear; clc;

% Time axis
T_total = 100;                 % total simulation time [s]
N       = 10001;               % time samples
t_s     = linspace(0, T_total, N);

% Breathing pattern
f_bpm   = 12 * ones(size(t_s));       % 12 breaths/min constant
VT_L    = 0.5 * ones(size(t_s));      % 0.5 L tidal volume
Vres_L  = 1.2 * ones(size(t_s));      % residual volume ~1.2 L (for later)
Vbase_L = 2.4 * ones(size(t_s));      % baseline ~FRC = 2.4 L

% Belly breathing fraction and inspiratory fraction
bellyFrac = 0.6 * ones(size(t_s));    % 60% belly breathing
inspFrac  = (1/3) * ones(size(t_s));  % I:E ≈ 1:2

lungParameters = struct( ...
    'f_bpm', f_bpm, ...
    'VT_L', VT_L, ...
    'Vres_L', Vres_L, ...
    'Vbase_L', Vbase_L, ...
    'bellyFrac', bellyFrac, ...
    'inspFrac', inspFrac, ...
    'lungSeparation_mm', 5 * ones(size(t_s)));

% Compute breathing motion + ellipsoid geometry
[leftEllipsoidParams, rightEllipsoidParams] = lung_ellipsoid_waveform(t_s, lungParameters);

% Visualization options
visOpts.frameStep = 4;        % skip frames for speed
visOpts.nTheta    = 200;

% Display breathing motion
visualize_lung_motion(t_s, leftEllipsoidParams, rightEllipsoidParams, visOpts);
