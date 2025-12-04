function [V_L, R_mm, H_mm, B_phase] = computeBreathingMotionEllipsoid( ...
            t_s, f_bpm, VT_L, Vres_L, Vbase_L, bellyFrac, inspFrac)
%COMPUTEBREATHINGMOTIONELLIPSOID  Lung volume + ellipsoid geometry vs time.
%
%   [V_L, R_mm, H_mm, B_phase] = computeBreathingMotionEllipsoid( ...
%           t_s, f_bpm, VT_L, Vres_L, Vbase_L, bellyFrac, inspFrac)
%
%   Inputs (1xN row vectors):
%       t_s       - time [s], strictly increasing
%       f_bpm     - breathing frequency [breaths/min]
%       VT_L      - tidal volume [L] (peak insp - base)
%       Vres_L    - residual volume [L] (not used in geometry, for later use)
%       Vbase_L   - baseline lung volume [L] (e.g. FRC) at each time
%       bellyFrac - belly-breathing fraction in [0,1] at each time
%       inspFrac  - inspiratory fraction of cycle in (0,1) at each time
%
%   Outputs:
%       V_L     - total lung volume [L] (both lungs combined)
%       R_mm    - effective semi-axis radius [mm] of each lung ellipsoid (LR/AP)
%       H_mm    - semi-axis height [mm] of each lung ellipsoid (SI)
%       B_phase - cumulative breathing phase [rad]
%
%   Model details:
%       - Breathing phase is built from instantaneous f_bpm(t).
%       - Within-cycle phase φ in [0,1) is mapped to a half-cosine insp/exp
%         waveform using local inspFrac(k).
%       - Volume: V_L(t) = Vbase_L(t) + VT_L(t) * g(φ, inspFrac).
%       - Geometry: two identical rotational ellipsoids:
%             V_total = (4/3) * pi * R^2 * H
%         Chest component (1-bellyFrac) changes R only, belly component
%         changes H only, applied sequentially so that V_total(t) is exact.

    arguments
        t_s       (1,:) double {mustBeReal, mustBeFinite}
        f_bpm     (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        VT_L      (1,:) double {mustBeReal, mustBeFinite}
        Vres_L    (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative} %#ok<INUSA>
        Vbase_L   (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        bellyFrac (1,:) double {mustBeReal, mustBeFinite}
        inspFrac  (1,:) double {mustBeReal, mustBeFinite}
    end

    N = numel(t_s);
    if any(diff(t_s) < 0)
        error('t_s must be strictly increasing.');
    end
    if any([numel(f_bpm), numel(VT_L), numel(Vres_L), ...
            numel(Vbase_L), numel(bellyFrac), numel(inspFrac)] ~= N)
        error('All input vectors must have the same length as t_s.');
    end
    if any(bellyFrac < 0 | bellyFrac > 1)
        error('bellyFrac must be in [0,1].');
    end

    % Clamp inspFrac to avoid exactly 0 or 1
    inspFrac = max(min(inspFrac, 0.99), 0.01);

    % ---------------------------------------------------------------------
    % 1) Breathing phase B_phase from instantaneous frequency f_bpm(t)
    % ---------------------------------------------------------------------
    f_Hz = f_bpm / 60;          % [Hz]
    B_phase = zeros(1, N);      % [rad]
    dt = diff(t_s);

    for k = 2:N
        B_phase(k) = B_phase(k-1) + 2*pi*f_Hz(k-1)*dt(k-1);
    end

    % Within-cycle phase φ in [0,1)
    phi_cycle = mod(B_phase, 2*pi) / (2*pi);   % dimensionless

    % ---------------------------------------------------------------------
    % 2) Dimensionless volume waveform g(φ, inspFrac) in [0,1]
    % ---------------------------------------------------------------------
    g = zeros(1, N);

    for k = 1:N
        alpha = inspFrac(k);        % local inspiratory fraction
        phi  = phi_cycle(k);

        if phi < alpha
            % Inspiration: 0 <= phi < alpha
            s = phi / alpha;        % [0,1]
            g(k) = 0.5*(1 - cos(pi*s));   % 0 -> 1
        else
            % Expiration: alpha <= phi < 1
            s = (phi - alpha) / (1 - alpha);  % [0,1]
            g(k) = 0.5*(1 + cos(pi*s));       % 1 -> 0
        end
    end

    % ---------------------------------------------------------------------
    % 3) Total lung volume V_L(t) [L]
    % ---------------------------------------------------------------------
    V_L = Vbase_L + VT_L .* g;

    % Convert to m^3
    V_m3 = V_L * 1e-3;

    % ---------------------------------------------------------------------
    % 4) Ellipsoid geometry: pair of identical lungs
    %
    %    Total volume for both lungs:
    %        V_tot = (4/3) * pi * R^2 * H
    %
    %    Baseline geometry: choose baseline radius R0_m, then compute H0_m
    %    from V_tot(1). This avoids needing H0_m as an explicit input.
    % ---------------------------------------------------------------------
    V0_m3 = V_m3(1);

    % Choose a baseline radius ~5 cm (can be adjusted)
    R0_m = 0.05;                               % 5 cm
    H0_m = (3 * V0_m3) / (4 * pi * R0_m^2);    % from V0 = (4/3)*pi*R0^2*H0

    % Total volume change relative to baseline
    dV_m3 = V_m3 - V0_m3;

    chestFrac = 1 - bellyFrac;
    dV_chest  = chestFrac .* dV_m3;
    dV_belly  = bellyFrac .* dV_m3;

    % ---------------------------------------------------------------------
    % 5) Apply chest component: change R only, keep H = H0
    %
    %    V1_tot = V0 + dV_chest = (4/3)*pi*R_ch^2*H0  =>  R_ch
    % ---------------------------------------------------------------------
    V1_m3 = V0_m3 + dV_chest;
    if any(V1_m3 <= 0)
        error('Intermediate volume became non-positive. Check inputs.');
    end

    R_ch_m = sqrt( (3 .* V1_m3) ./ (4 * pi * H0_m) );

    % ---------------------------------------------------------------------
    % 6) Apply belly component: change H only, keep R = R_ch
    %
    %    V2_tot = V1 + dV_belly = (4/3)*pi*R_ch^2*H  =>  H
    % ---------------------------------------------------------------------
    V2_m3 = V1_m3 + dV_belly;
    if any(V2_m3 <= 0)
        error('Final volume became non-positive. Check inputs.');
    end

    H_m = (3 .* V2_m3) ./ (4 * pi .* R_ch_m.^2);
    R_m = R_ch_m;

    % Convert to millimeters for output
    R_mm = 1000 .* R_m;
    H_mm = 1000 .* H_m;

    % Optional sanity check:
    % V_check = (4/3)*pi.*R_m.^2.*H_m;
    % maxRelErr = max(abs(V_check - V_m3)./max(V_m3, eps));
end
