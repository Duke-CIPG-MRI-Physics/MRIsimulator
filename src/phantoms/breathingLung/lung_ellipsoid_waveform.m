function [R_mm, H_mm] = lung_ellipsoid_waveform(t_s, lungParameters)
%LUNG_ELLIPSOID_WAVEFORM  Lung ellipsoid geometry vs time.
%
%   [leftEllipsoidParams, rightEllipsoidParams] = lung_ellipsoid_waveform(t_s, lungParameters)
%
%   Inputs:
%       t_s     - time [s], strictly increasing 1xN vector
%       lungParameters - struct with 1xN vectors for fields:
%           f_bpm            - breathing frequency [breaths/min]
%           VT_L             - tidal volume [L] (peak insp - base)
%           Vres_L           - residual volume [L] (kept for compatibility)
%           Vbase_L          - baseline lung volume [L] (e.g. FRC)
%           bellyFrac        - belly-breathing fraction in [0,1]
%           inspFrac         - inspiratory fraction of cycle in (0,1)
%           lungSeparation_mm- spacing added to the lung radius [mm]
%
%   Outputs:
%       leftEllipsoidParams, rightEllipsoidParams - struct with fields:
%           .a_mm, .b_mm, .c_mm  - ellipsoid semi-axes [mm]
%           .pose               - pose struct with center and orientation
%         Additional convenience fields include the underlying waveform
%         values (.R_mm, .H_mm, .lungPosition_mm).
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
        t_s             (1,:) double {mustBeReal, mustBeFinite}
        lungParameters  struct
    end

    requiredFields = {'f_bpm', 'VT_L', 'Vres_L', 'Vbase_L', ...
        'bellyFrac', 'inspFrac'};
    for idxField = 1:numel(requiredFields)
        if ~isfield(lungParameters, requiredFields{idxField})
            error('lung_ellipsoid_waveform:MissingField', ...
                'lungParameters.%s is required.', requiredFields{idxField});
        end
    end

    f_bpm = lungParameters.f_bpm;
    VT_L  = lungParameters.VT_L;
    Vres_L = lungParameters.Vres_L; %#ok<NASGU>
    Vbase_L = lungParameters.Vbase_L;
    bellyFrac = lungParameters.bellyFrac;
    inspFrac  = lungParameters.inspFrac;

    validateattributes(f_bpm, {'numeric'}, {'real', 'finite', 'nonnegative'});
    validateattributes(VT_L, {'numeric'}, {'real', 'finite'});
    validateattributes(Vres_L, {'numeric'}, {'real', 'finite', 'nonnegative'});
    validateattributes(Vbase_L, {'numeric'}, {'real', 'finite', 'positive'});
    validateattributes(bellyFrac, {'numeric'}, {'real', 'finite'});
    validateattributes(inspFrac, {'numeric'}, {'real', 'finite'});
    
    N = numel(t_s);
    if any(diff(t_s) < 0)
        error('t_s must be strictly increasing.');
    end

    paramNames = {'f_bpm', 'VT_L', 'Vres_L', 'Vbase_L', ...
        'bellyFrac', 'inspFrac', 'lungSeparation_mm'};
    paramValues = {f_bpm, VT_L, Vres_L, Vbase_L, bellyFrac, inspFrac};

    nonScalarSize = [];
    for idxParam = 1:numel(paramValues)
        thisValue = paramValues{idxParam};
        if isscalar(thisValue)
            continue
        end

        thisSize = size(thisValue);
        if isempty(nonScalarSize)
            if numel(thisValue) ~= N
                error('lung_ellipsoid_waveform:LengthMismatch', ...
                    ['Non-scalar lungParameters.%s must have %d elements ' ...
                    'to align with t_s.'], paramNames{idxParam}, N);
            end
            nonScalarSize = thisSize;
        elseif ~isequal(thisSize, nonScalarSize)
            error('lung_ellipsoid_waveform:SizeMismatch', ...
                ['All non-scalar lungParameters must have the same size. ' ...
                '%s has size %s; expected %s.'], paramNames{idxParam}, ...
                mat2str(thisSize), mat2str(nonScalarSize));
        end
    end

    if isempty(nonScalarSize)
        nonScalarSize = size(t_s);
    end

    if any(bellyFrac < 0 | bellyFrac > 1)
        error('bellyFrac must be in [0,1].');
    end

    % Clamp inspFrac to avoid exactly 0 or 1
    inspFrac = max(min(inspFrac, 0.99), 0.01);

    %% Expand scalars for vectorized evaluation
    if isscalar(f_bpm)
        f_bpm = repmat(f_bpm, nonScalarSize);
    end
    if isscalar(VT_L)
        VT_L = repmat(VT_L, nonScalarSize);
    end
    if isscalar(Vbase_L)
        Vbase_L = repmat(Vbase_L, nonScalarSize);
    end
    if isscalar(bellyFrac)
        bellyFrac = repmat(bellyFrac, nonScalarSize);
    end
    if isscalar(inspFrac)
        inspFrac = repmat(inspFrac, nonScalarSize);
    end

    % ---------------------------------------------------------------------
    % 1) Breathing phase B_phase from instantaneous frequency f_bpm(t)
    % ---------------------------------------------------------------------
    waveformSize = nonScalarSize;

    f_Hz = f_bpm / 60;          % [Hz]
    dt = diff(t_s(:));
    f_Hz_vec = f_Hz(:);
    B_phase_vec = [0; cumsum(2*pi*f_Hz_vec(1:end-1) .* dt)];
    B_phase = reshape(B_phase_vec, waveformSize);      % [rad]

    % Within-cycle phase φ in [0,1)
    phi_cycle = mod(B_phase, 2*pi) / (2*pi);   % dimensionless

    % ---------------------------------------------------------------------
    % 2) Dimensionless volume waveform g(φ, inspFrac) in [0,1]
    % ---------------------------------------------------------------------
    phi_vec = phi_cycle(:);
    alpha_vec = inspFrac(:);
    g_vec = zeros(numel(phi_vec), 1);
    isInspiration = phi_vec < alpha_vec;
    s = zeros(numel(phi_vec), 1);
    s(isInspiration) = phi_vec(isInspiration) ./ alpha_vec(isInspiration);
    g_vec(isInspiration) = 0.5 * (1 - cos(pi * s(isInspiration)));
    s(~isInspiration) = (phi_vec(~isInspiration) - alpha_vec(~isInspiration)) ...
        ./ (1 - alpha_vec(~isInspiration));
    g_vec(~isInspiration) = 0.5 * (1 + cos(pi * s(~isInspiration)));
    g = reshape(g_vec, waveformSize);

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
