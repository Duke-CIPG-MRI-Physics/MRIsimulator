% demo_cardiac_waveform.m
% Demo of time-varying LV ellipsoid geometry + volume + strain visualization
% The core waveform model now evaluates long time vectors in 15,625-sample
% chunks (~0.001 Gb of temporaries) to limit memory use.
%
% Requires:
%   - cardiac_ellipsoid_waveform.m
%   - visualize_cardiac_ellipsoid.m

clear; clc; close all;

%% 1) Time axis and basic hemodynamics
T_total_s = 6;                    % total simulation time [s]
N         = 2001;                 % number of time samples
t_s       = linspace(0, T_total_s, N);

% Heart rate: 70 bpm constant (you can make this vary over time if desired)
HR_bpm = 70 * ones(1, N);

% LV volumes: EDV/ESV constant over time for this demo
EDV_ml = 120 * ones(1, N);        % end-diastolic volume [mL]
ESV_ml = 50  * ones(1, N);        % end-systolic volume [mL]

SV_ml = EDV_ml(1) - ESV_ml(1);    % stroke volume [mL]
EF    = SV_ml / EDV_ml(1);        % ejection fraction

fprintf('Demo EF = %.1f %% (SV = %.1f mL, EDV = %.1f mL)\n', 100*EF, SV_ml, EDV_ml(1));

%% 2) Waveform + geometry options
opts = struct();
opts.systFrac  = 0.35;    % systolic fraction of cycle
opts.q_ED      = 2.5;     % ED aspect ratio a0/b0 (long/short)
opts.GLS_peak  = -0.20;   % peak longitudinal strain (GLS)
opts.GCS_peak  = -0.25;   % peak circumferential strain (GCS)

%% 3) Run ellipsoid waveform model
[V_ml, a_mm, b_mm, c_mm, phase, eps_L, eps_C] = ...
    cardiac_ellipsoid_waveform(t_s, HR_bpm, EDV_ml, ESV_ml, opts);

% Optional: equivalent sphere radius (sanity check)
V_m3   = V_ml * 1e-6;
R_eq_m = (3*V_m3 ./ (4*pi)).^(1/3); %#ok<NASGU>

%% 4) Visualization options and animation
visOpts = struct();
visOpts.frameStep = 4;     % subsample frames to speed up animation (~N/4 frames)
visOpts.nTheta    = 200;   % angular resolution for ellipse/circle

visualize_cardiac_ellipsoid(t_s, V_ml, a_mm, b_mm, c_mm, eps_L, eps_C, visOpts);
