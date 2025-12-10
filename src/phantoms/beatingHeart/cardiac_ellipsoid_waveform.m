function [V_ml, a_mm, b_mm, c_mm, phase, eps_L, eps_C] = ...
    cardiac_ellipsoid_waveform(t_s, HR_bpm, EDV_ml, ESV_ml, opts)
%CARDIAC_ELLIPSOID_WAVEFORM  LV volume + ellipsoid semi-axes vs time.
%
%   [V_ml, a_mm, b_mm, c_mm, phase, eps_L, eps_C] = ...
%       cardiac_ellipsoid_waveform(t_s, HR_bpm, EDV_ml, ESV_ml, opts)
%
%   Inputs (1xN vectors or scalars used directly):
%       t_s     - time [s], strictly increasing
%       HR_bpm  - heart rate [beats/min] at each time
%       EDV_ml  - LV end-diastolic volume [mL] at each time (max volume)
%       ESV_ml  - LV end-systolic volume [mL] at each time (min volume)
%
%   opts (struct, optional fields):
%       .systFrac  - systolic fraction of R-R interval (0<beta<1), default 0.35
%       .q_ED      - ED aspect ratio a0/b0 (long/short), default 2.5
%       .GLS_peak  - peak global longitudinal strain (negative), default -0.20
%       .GCS_peak  - peak global circumferential strain (negative), default -0.25
%
%   Outputs:
%       V_ml   - instantaneous LV volume [mL] vs time
%       a_mm   - ellipsoid semi-axis (radius) along LV long-axis [mm] vs time
%       b_mm   - ellipsoid short semi-axis (radius) [mm] vs time
%       c_mm   - ellipsoid short semi-axis (radius, =b_mm) [mm] vs time
%       phase  - cumulative cardiac phase [rad] (0..∞), built from HR(t)
%       eps_L  - longitudinal strain vs time (relative to ED)
%       eps_C  - circumferential strain vs time (relative to ED)
%
%   Model:
%       1) Build cumulative phase from instantaneous HR(t).
%       2) Build LV volume V(t) with a piecewise half-cosine EDV->ESV->EDV.
%       3) Determine ED ellipsoid geometry (a0,b0,c0) from EDV and q_ED.
%       4) Build strain-driven template axes (ahat,bhat,chat) using GLS/GCS.
%       5) At each t, rescale template axes isotropically so that
%              V(t) = (4/3)*pi*a(t)*b(t)*c(t)
%          exactly, preserving your volume waveform and EF.

    arguments
        t_s    (1,:) double {mustBeReal, mustBeFinite}
        HR_bpm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        EDV_ml (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        ESV_ml (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        opts struct = struct()
    end

    % ---------------------------------------------------------------------
    % Defaults for opts
    % ---------------------------------------------------------------------
    if ~isfield(opts, 'systFrac') || isempty(opts.systFrac)
        opts.systFrac = 0.35;
    end
    if ~isfield(opts, 'q_ED') || isempty(opts.q_ED)
        opts.q_ED = 2.5;
    end
    if ~isfield(opts, 'GLS_peak') || isempty(opts.GLS_peak)
        opts.GLS_peak = -0.20;
    end
    if ~isfield(opts, 'GCS_peak') || isempty(opts.GCS_peak)
        opts.GCS_peak = -0.25;
    end

    beta     = opts.systFrac;
    q        = opts.q_ED;
    GLS_peak = opts.GLS_peak;
    GCS_peak = opts.GCS_peak;

    if beta <= 0 || beta >= 1
        error('opts.systFrac must be in (0,1).');
    end

    % ---------------------------------------------------------------------
    % Basic checks
    % ---------------------------------------------------------------------
    t_s = t_s(:)';
    N = numel(t_s);
    if any(diff(t_s) < 0)
        error('t_s must be strictly increasing.');
    end

    HR_bpm = validateLengthOrScalar(HR_bpm, N, 'HR_bpm');
    EDV_ml = validateLengthOrScalar(EDV_ml, N, 'EDV_ml');
    ESV_ml = validateLengthOrScalar(ESV_ml, N, 'ESV_ml');

    if any(ESV_ml >= EDV_ml)
        error('ESV_ml must be < EDV_ml at all time points.');
    end

    % ---------------------------------------------------------------------
    % 1) Build cumulative phase from HR(t)
    % ---------------------------------------------------------------------
    f_Hz = HR_bpm / 60;     % [Hz]
    dt = diff(t_s);         % N-1 intervals for N samples

    % Use HR samples per interval (N-1). For scalar HR, MATLAB broadcasts the
    % multiplication; for vectors we trim the unused last element to match dt.
    f_interval = f_Hz;
    if ~isscalar(f_interval)
        f_interval = f_interval(1:end-1);
    end

    incrementalPhase = 2*pi .* f_interval .* dt;
    phase = [0 cumsum(incrementalPhase)];
    % for k = 2:N
    %     phase(k) = phase(k-1) + 2*pi*f_Hz(k-1)*dt(k-1);
    % end

    % Within-cycle phase φ in [0,1)
    phi_cycle = mod(phase, 2*pi) / (2*pi);   % dimensionless

    % ---------------------------------------------------------------------
    % 2) Volume waveform via piecewise half-cosines
    % ---------------------------------------------------------------------
    V_ml = zeros(1, N);
    SV_ml = EDV_ml - ESV_ml;

    systMask = (phi_cycle < beta);
    diasMask = ~systMask;

    % Systole: EDV -> ESV (volume falls)
    s_syst = phi_cycle(systMask) / beta;       % [0,1]
    h_syst = 0.5*(1 + cos(pi*s_syst));         % 1 -> 0

    % Diastole: ESV -> EDV (volume rises)
    s_dias = (phi_cycle(diasMask) - beta) / (1 - beta);  % [0,1]
    h_dias = 0.5*(1 - cos(pi*s_dias));                   % 0 -> 1

    waveform = zeros(1, N);
    waveform(systMask) = h_syst;
    waveform(diasMask) = h_dias;

    V_ml = ESV_ml + SV_ml .* waveform;

    % Convert volume to m^3
    V_m3 = V_ml * 1e-6;

    % ---------------------------------------------------------------------
    % 3) Baseline ED ellipsoid from EDV and aspect ratio q = a0/b0
    % ---------------------------------------------------------------------
    [~, idxED] = max(V_m3);       % index of ED
    V_ED = V_m3(idxED);           % ED volume [m^3]

    % Prolate spheroid: a0 = q * b0, c0 = b0
    % V_ED = (4/3)*pi*a0*b0^2 = (4/3)*pi*q*b0^3  =>  solve for b0
    b0 = ((3*V_ED) / (4*pi*q))^(1/3);
    a0 = q * b0;
    % c0 = b0;  % implied

    % ---------------------------------------------------------------------
    % 4) Strain-driven template axes (relative to ED)
    % ---------------------------------------------------------------------
    eps_L = zeros(1, N);   % longitudinal strain
    eps_C = zeros(1, N);   % circumferential strain

    % Reuse systolic/diastolic masks
    systMask = (phi_cycle < beta);
    diasMask = ~systMask;

    % Systolic part: 0→beta
    s_syst = phi_cycle(systMask) / beta;            % [0,1]
    f_syst = 0.5*(1 - cos(pi*s_syst));              % 0 -> 1
    eps_L(systMask) = GLS_peak * f_syst;
    eps_C(systMask) = GCS_peak * f_syst;

    % Diastolic part: beta→1
    s_dias = (phi_cycle(diasMask) - beta) / (1 - beta);  % [0,1]
    f_dias = 0.5*(1 + cos(pi*s_dias));                   % 1 -> 0
    eps_L(diasMask) = GLS_peak * f_dias;
    eps_C(diasMask) = GCS_peak * f_dias;

    % Template semi-axes (before volume correction)
    a_hat = a0 * (1 + eps_L);
    b_hat = b0 * (1 + eps_C);
    c_hat = b_hat;                % rotationally symmetric short axes

    % Template volume
    V_hat_m3 = (4/3)*pi .* a_hat .* b_hat.^2;

    % ---------------------------------------------------------------------
    % 5) Isotropic scale factor to enforce exact V(t)
    % ---------------------------------------------------------------------
    if any(V_hat_m3 <= 0)
        error('Template ellipsoid volume became non-positive; check strain settings.');
    end
    lambda = (V_m3 ./ V_hat_m3).^(1/3);

    % Final semi-axes (meters)
    a_m = lambda .* a_hat;
    b_m = lambda .* b_hat;
    c_m = lambda .* c_hat;

    % Convert to millimeters for output
    a_mm = 1000 .* a_m;
    b_mm = 1000 .* b_m;
    c_mm = 1000 .* c_m;

    % Optional internal check (commented out)
    % V_check = (4/3)*pi .* a_m .* b_m.^2;
    % maxRelErr = max(abs(V_check - V_m3)./max(V_m3, eps));
    % fprintf('Max relative volume error = %.3g\n', maxRelErr);

end

% -------------------------------------------------------------------------
function vec = validateLengthOrScalar(vec, N, name)
%VALIDATELENGTHORSCALAR  Ensure input is a row vector with numel 1 or N.

    vec = vec(:)';
    if ~(isscalar(vec) || numel(vec) == N)
        error('%s must be scalar or have the same number of elements as t_s.', name);
    end
end
