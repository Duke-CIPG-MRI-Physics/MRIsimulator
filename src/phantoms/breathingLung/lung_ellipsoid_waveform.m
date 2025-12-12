function [leftEllipsoidParams, rightEllipsoidParams] = lung_ellipsoid_waveform(t_s, lungParameters)
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
        'bellyFrac', 'inspFrac', 'lungSeparation_mm'};
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
    lungSeparation_mm = lungParameters.lungSeparation_mm;

    validateattributes(f_bpm, {'numeric'}, {'real', 'finite', 'nonnegative'});
    validateattributes(VT_L, {'numeric'}, {'real', 'finite'});
    validateattributes(Vres_L, {'numeric'}, {'real', 'finite', 'nonnegative'});
    validateattributes(Vbase_L, {'numeric'}, {'real', 'finite', 'positive'});
    validateattributes(bellyFrac, {'numeric'}, {'real', 'finite'});
    validateattributes(inspFrac, {'numeric'}, {'real', 'finite'});
    validateattributes(lungSeparation_mm, {'numeric'}, {'real', 'finite', 'nonnegative'});

    N = numel(t_s);
    if any(diff(t_s) < 0)
        error('t_s must be strictly increasing.');
    end

    paramNames = {'f_bpm', 'VT_L', 'Vres_L', 'Vbase_L', ...
        'bellyFrac', 'inspFrac', 'lungSeparation_mm'};
    paramValues = {f_bpm, VT_L, Vres_L, Vbase_L, bellyFrac, inspFrac, ...
        lungSeparation_mm};

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

    % ---------------------------------------------------------------------
    % 1) Breathing phase B_phase from instantaneous frequency f_bpm(t)
    % ---------------------------------------------------------------------
    waveformSize = nonScalarSize;

    f_Hz = f_bpm / 60;          % [Hz]
    B_phase = zeros(waveformSize);      % [rad]
    dt = diff(t_s);

    for k = 2:N
        f_Hz_km1 = valueAtIndex(f_Hz, k-1);
        B_phase(k) = B_phase(k-1) + 2*pi*f_Hz_km1*dt(k-1);
    end

    % Within-cycle phase φ in [0,1)
    phi_cycle = mod(B_phase, 2*pi) / (2*pi);   % dimensionless

    % ---------------------------------------------------------------------
    % 2) Dimensionless volume waveform g(φ, inspFrac) in [0,1]
    % ---------------------------------------------------------------------
    g = zeros(waveformSize);

    for k = 1:N
        alpha = valueAtIndex(inspFrac, k);        % local inspiratory fraction
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

    lungPosition_mm = R_mm + lungSeparation_mm;

    lungParams = struct('a_mm', R_mm, ...
        'b_mm', R_mm, ...
        'c_mm', H_mm, ...
        'R_mm', R_mm, ...
        'H_mm', H_mm, ...
        'lungPosition_mm', lungPosition_mm);

    rightCenter_mm = [lungPosition_mm(:), ...
        zeros(numel(lungPosition_mm), 1), ...
        zeros(numel(lungPosition_mm), 1)];
    leftCenter_mm = [-lungPosition_mm(:), ...
        zeros(numel(lungPosition_mm), 1), ...
        zeros(numel(lungPosition_mm), 1)];

    rightEllipsoidParams = lungParams;
    rightEllipsoidParams.pose = struct('center', struct('x_mm', rightCenter_mm(:,1), ...
        'y_mm', rightCenter_mm(:,2), ...
        'z_mm', rightCenter_mm(:,3)), ...
        'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);
    leftEllipsoidParams = lungParams;
    leftEllipsoidParams.pose = struct('center', struct('x_mm', leftCenter_mm(:,1), ...
        'y_mm', leftCenter_mm(:,2), ...
        'z_mm', leftCenter_mm(:,3)), ...
        'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);

    % Optional sanity check:
    % V_check = (4/3)*pi.*R_m.^2.*H_m;
    % maxRelErr = max(abs(V_check - V_m3)./max(V_m3, eps));
end

% -------------------------------------------------------------------------
function value = valueAtIndex(param, idx)
%VALUEATINDEX Return a scalar parameter value, tolerating scalar inputs.
    if isscalar(param)
        value = param;
    else
        value = param(idx);
    end
end
