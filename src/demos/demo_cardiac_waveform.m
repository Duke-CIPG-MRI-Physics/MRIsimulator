% demo_cardiac_waveform.m
% Demo of time-varying LV ellipsoid geometry using the simplified interface
%   of cardiac_ellipsoid_waveform (returns a, b, and phase only).
% The core waveform model now evaluates long time vectors in 15,625-sample
% chunks (~0.001 Gb of temporaries) to limit memory use.
%
% Requires:
%   - cardiac_ellipsoid_waveform.m

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

%% 3) Run ellipsoid waveform model (new interface: a, b, phase)
[a_mm, b_mm, phase] = cardiac_ellipsoid_waveform(t_s, HR_bpm, EDV_ml, ESV_ml, opts); %#ok<NASGU>

%% 4) Visualization: ellipse outline + semi-axes vs time
frameStep = 4;             % subsample frames to speed up animation (~N/4 frames)
nTheta    = 200;           % angular resolution for ellipse

a_cm = a_mm / 10;          % convert to cm for plotting
b_cm = b_mm / 10;

theta = linspace(0, 2*pi, nTheta);

% Axis limits
maxA = max(a_cm);
maxB = max(b_cm);
pad = 1.2;

figure('Name','Cardiac ellipsoid (simplified)','Color','w');
tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% (1) Coronal ellipse outline (x = b, y = a)
axEll = nexttile(tlo, 1);
hold(axEll, 'on');
x0 = b_cm(1) * cos(theta);
y0 = a_cm(1) * sin(theta);
hEll = plot(axEll, x0, y0, 'b', 'LineWidth', 2);
axis(axEll, 'equal');
xlim(axEll, pad * [-maxB, maxB]);
ylim(axEll, pad * [-maxA, maxA]);
xlabel(axEll, 'Left–Right [cm]');
ylabel(axEll, 'Long-axis [cm]');
title(axEll, 'Ellipsoid cross-section');
grid(axEll, 'on');

% (2) Semi-axes vs time with moving cursor
axAxes = nexttile(tlo, 2);
hold(axAxes, 'on');
plot(axAxes, t_s, a_cm, 'LineWidth', 1.5, 'DisplayName', 'a (long)');
plot(axAxes, t_s, b_cm, 'LineWidth', 1.5, 'DisplayName', 'b (short)');
hPointA = plot(axAxes, t_s(1), a_cm(1), 'ok', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 6, 'DisplayName', 'Cursor a');
hPointB = plot(axAxes, t_s(1), b_cm(1), 'or', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 6, 'DisplayName', 'Cursor b');
xlabel(axAxes, 'Time [s]');
ylabel(axAxes, 'Semi-axis [cm]');
title(axAxes, 'Semi-axes vs time');
legend(axAxes, 'Location', 'best');
grid(axAxes, 'on');

% Animation loop
for k = 1:frameStep:numel(t_s)
    ak = a_cm(k);
    bk = b_cm(k);

    xk = bk * cos(theta);
    yk = ak * sin(theta);
    set(hEll, 'XData', xk, 'YData', yk);

    set(hPointA, 'XData', t_s(k), 'YData', ak);
    set(hPointB, 'XData', t_s(k), 'YData', bk);

    drawnow;
end