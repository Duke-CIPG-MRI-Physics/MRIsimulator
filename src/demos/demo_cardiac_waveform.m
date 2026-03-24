% demo_cardiac_waveform.m
% Demo of the time-varying LV ellipsoid waveform using the current
% cardiac_ellipsoid_waveform interface.
%
% The cardiac waveform helper returns ellipsoid parameters in the phantom
% axis convention, where the long-axis radius is stored in c_mm. This demo
% derives the LV volume and axis-based strains from that geometry, then
% uses visualize_cardiac_ellipsoid to show the shape and waveform over
% time.

clear; clc; close all;

%% 1) Time axis and waveform options
T_total_s = 6;
N = 2001;
t_s = linspace(0, T_total_s, N);

cardiacOpts = struct();
cardiacOpts.systFrac = 0.35;
cardiacOpts.q_ED = 2.5;
cardiacOpts.GLS_peak = -0.20;
cardiacOpts.GCS_peak = -0.25;
cardiacOpts.HR_bpm = 70;
cardiacOpts.EDV_ml = 120;
cardiacOpts.ESV_ml = 50;

%% 2) Run the waveform model and derive plotted quantities
ellipsoidParams = cardiac_ellipsoid_waveform(t_s, cardiacOpts);

shortAxisRadius_mm = ellipsoidParams.b_mm(1, :);
longAxisRadius_mm = ellipsoidParams.c_mm(1, :);
volume_ml = (4/3) * pi .* ellipsoidParams.a_mm(1, :) .* ...
    ellipsoidParams.b_mm(1, :) .* ellipsoidParams.c_mm(1, :) / 1000;

endDiastolicLongAxis_mm = max(longAxisRadius_mm);
endDiastolicShortAxis_mm = max(shortAxisRadius_mm);
eps_L = longAxisRadius_mm ./ endDiastolicLongAxis_mm - 1;
eps_C = shortAxisRadius_mm ./ endDiastolicShortAxis_mm - 1;

strokeVolume_ml = max(volume_ml) - min(volume_ml);
ejectionFraction = strokeVolume_ml / max(volume_ml);

fprintf('Demo EF = %.1f %% (SV = %.1f mL, EDV = %.1f mL)\n', ...
    100 * ejectionFraction, strokeVolume_ml, max(volume_ml));

%% 3) Visualize waveform and ellipsoid geometry over time
visOpts = struct();
visOpts.frameStep = 4;
visOpts.nTheta = 200;

visualize_cardiac_ellipsoid(t_s, volume_ml, longAxisRadius_mm, ...
    shortAxisRadius_mm, shortAxisRadius_mm, eps_L, eps_C, visOpts);
