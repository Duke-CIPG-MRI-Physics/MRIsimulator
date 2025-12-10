function [a_mm, b_mm, phase] = ...
    cardiac_ellipsoid_waveform(t_s, HR_bpm, EDV_ml, ESV_ml, opts)
%CARDIAC_ELLIPSOID_WAVEFORM  LV ellipsoid semi-axes driven by volume + strain.
%
%   [a_mm, b_mm, phase] = cardiac_ellipsoid_waveform(t_s, HR_bpm, EDV_ml,
%       ESV_ml, opts)
%
%   Inputs (1xN vectors or scalars broadcasted to length N):
%       t_s     - time [s], strictly increasing (length N)
%       HR_bpm  - heart rate [beats/min] at each time sample
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
%       a_mm   - ellipsoid semi-axis (radius) along LV long-axis [mm] vs time
%       b_mm   - ellipsoid short semi-axis (radius) [mm] vs time
%       phase  - cumulative cardiac phase [cycles], built from HR(t)
%
%   Model:
%       1) Build cumulative phase from instantaneous HR(t).
%       2) Build LV volume V(t) with a piecewise half-cosine EDV->ESV->EDV.
%       3) Determine ED ellipsoid geometry (a0,b0,c0) from EDV and q_ED.
%       4) Build strain-driven template axes (ahat,bhat,chat) using GLS/GCS.
%       5) At each t, rescale template axes isotropically so that
%              V(t) = (4/3)*pi*a(t)*b(t)*c(t)
%          exactly, preserving your volume waveform and EF.
%
%   Memory usage:
%       After computing phase, volume and strain calculations are processed in
%       chunks of up to 15,625 samples (~0.001 Gb of doubles) to avoid large
%       temporary matrices. Shorter signals are evaluated directly without
%       chunking so the chunk calculations minimally affect memory use.

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
    numSamples = numel(t_s);
    dt = diff(t_s(:)');
    if any(dt <= 0)
        error('t_s must be strictly increasing.');
    end

    HR_bpm = validateLengthOrScalar(HR_bpm, numSamples, 'HR_bpm', numSamples - 1);
    EDV_ml = validateLengthOrScalar(EDV_ml, numSamples, 'EDV_ml');
    ESV_ml = validateLengthOrScalar(ESV_ml, numSamples, 'ESV_ml');

    if any(ESV_ml >= EDV_ml)
        error('ESV_ml must be < EDV_ml at all time points.');
    end

    % ---------------------------------------------------------------------
    % 1) Build cumulative phase from HR(t)
    % ---------------------------------------------------------------------
    f_Hz = HR_bpm / 60;     % [Hz]

    % Use HR samples per interval (N-1). For scalar HR, MATLAB broadcasts the
    % multiplication; for vectors we trim the unused last element to match dt.
    f_interval = buildIntervalFrequencies(f_Hz, numel(dt), numSamples);

    phase = mod([0 cumsum(2*pi .* f_interval .* dt)], 2*pi) / (2*pi);
    clear dt f_interval;
    
    % ---------------------------------------------------------------------
    % 2) Volume waveform via piecewise half-cosines
    % ---------------------------------------------------------------------
    % Chunk calculations to ~0.001 Gb of temporary data so memory impact is
    % minimal.
    chunkSize = 15625;
    V_ED = -Inf;

    for idxStart = 1:chunkSize:numSamples
        idxEnd = min(idxStart + chunkSize - 1, numSamples);
        idx = idxStart:idxEnd;

        waveformChunk = zeros(1, numel(idx));
        eps_L_chunk = zeros(1, numel(idx));
        eps_C_chunk = zeros(1, numel(idx));

        [waveformChunk, eps_L_chunk, eps_C_chunk] = buildWaveformAndStrain( ...
            phase(idx), beta, GLS_peak, GCS_peak, waveformChunk, ...
            eps_L_chunk, eps_C_chunk);

        V_chunk = computeVolumeChunk(EDV_ml, ESV_ml, idx, waveformChunk);

        V_chunk_max = max(V_chunk);
        if V_chunk_max > V_ED
            V_ED = V_chunk_max;
        end
    end

    % ---------------------------------------------------------------------
    % 3) Baseline ED ellipsoid from EDV and aspect ratio q = a0/b0
    % ---------------------------------------------------------------------

    % Prolate spheroid: a0 = q * b0, c0 = b0
    % V_ED = (4/3)*pi*a0*b0^2 = (4/3)*pi*q*b0^3  =>  solve for b0
    b0 = ((3*V_ED) / (4*pi*q))^(1/3);
    a0 = q * b0;
    % c0 = b0;  % implied

    % ---------------------------------------------------------------------
    % 4-5) Strain-driven axes + isotropic scale factor (chunked)
    % ---------------------------------------------------------------------
    a_mm = zeros(1, numSamples);
    b_mm = zeros(1, numSamples);

    for idxStart = 1:chunkSize:numSamples
        idxEnd = min(idxStart + chunkSize - 1, numSamples);
        idx = idxStart:idxEnd;

        waveformChunk = zeros(1, numel(idx));
        eps_L_chunk = zeros(1, numel(idx));
        eps_C_chunk = zeros(1, numel(idx));

        [waveformChunk, eps_L_chunk, eps_C_chunk] = buildWaveformAndStrain( ...
            phase(idx), beta, GLS_peak, GCS_peak, waveformChunk, ...
            eps_L_chunk, eps_C_chunk);

        V_chunk = computeVolumeChunk(EDV_ml, ESV_ml, idx, waveformChunk);

        a_hat_chunk = a0 * (1 + eps_L_chunk);
        b_hat_chunk = b0 * (1 + eps_C_chunk);
        V_hat_chunk = (4/3)*pi .* a_hat_chunk .* b_hat_chunk.^2;

        if any(V_hat_chunk <= 0)
            error('Template ellipsoid volume became non-positive; check strain settings.');
        end

        lambda_chunk = (V_chunk ./ V_hat_chunk).^(1/3);

        a_mm(idx) = 1000 .* lambda_chunk .* a_hat_chunk;
        b_mm(idx) = 1000 .* lambda_chunk .* b_hat_chunk;
    end
end

% -------------------------------------------------------------------------
function vec = validateLengthOrScalar(vec, primaryLength, name, altLength)
%VALIDATELENGTHORSCALAR  Ensure input is a row vector with accepted lengths.

    if nargin < 4
        altLength = [];
    end

    vec = vec(:)';
    validLengths = [primaryLength, altLength];
    validLengths = validLengths(validLengths > 0);
    if ~(isscalar(vec) || any(numel(vec) == validLengths))
        if isempty(altLength)
            error('%s must be scalar or have %d elements.', name, primaryLength);
        else
            error('%s must be scalar or have %d elements (or %d if specified per interval).', ...
                name, primaryLength, primaryLength - 1);
        end
    end
end

% -------------------------------------------------------------------------
function f_interval = buildIntervalFrequencies(f_Hz, numIntervals, numSamples)
%BUILDINTERVALFREQUENCIES  Align heart rate waveform with interval durations.

    if isscalar(f_Hz)
        f_interval = f_Hz * ones(1, numIntervals);
        return;
    end

    numSamplesProvided = numel(f_Hz);
    if numSamplesProvided == numIntervals
        f_interval = f_Hz;
    elseif numSamplesProvided == numSamples
        f_interval = f_Hz(1:end-1);
    else
        error(['HR_bpm must be scalar or have the same length as t_s or diff(t_s). ', ...
            'Received %d samples.'], numSamplesProvided);
    end
end

% -------------------------------------------------------------------------
function V_chunk = computeVolumeChunk(EDV_ml, ESV_ml, idx, waveformChunk)
%COMPUTEVOLUMECHUNK  Volume waveform in cubic meters for a chunk of indices.

    EDV_slice = extractParamSlice(EDV_ml, idx);
    ESV_slice = extractParamSlice(ESV_ml, idx);
    SV_ml = EDV_slice - ESV_slice;
    V_chunk = (ESV_slice + SV_ml .* waveformChunk) * 1e-6;
end

% -------------------------------------------------------------------------
function [waveformChunk, eps_L_chunk, eps_C_chunk] = buildWaveformAndStrain( ...
        phaseChunk, beta, GLS_peak, GCS_peak, waveformChunk, eps_L_chunk, ...
        eps_C_chunk)
%BUILDWAVEFORMANDSTRAIN  Half-cosine waveform + strains for a phase chunk.

    systMask = (phaseChunk < beta);
    diasMask = ~systMask;

    if any(systMask)
        f_syst = 0.5*(1 - cos(pi*phaseChunk(systMask) / beta));              % 0 -> 1
        waveformChunk(systMask) = 0.5*(1 + cos(pi*phaseChunk(systMask) / beta));
        eps_L_chunk(systMask) = GLS_peak * f_syst;
        eps_C_chunk(systMask) = GCS_peak * f_syst;
    end

    if any(diasMask)
        shiftedPhase = phaseChunk(diasMask) - beta;
        f_dias = 0.5*(1 + cos(pi * shiftedPhase / (1 - beta)));                   % 1 -> 0
        waveformChunk(diasMask) = 0.5*(1 - cos(pi * shiftedPhase / (1 - beta)));
        eps_L_chunk(diasMask) = GLS_peak * f_dias;
        eps_C_chunk(diasMask) = GCS_peak * f_dias;
    end
end

% -------------------------------------------------------------------------
function slice = extractParamSlice(param, idx)
%EXTRACTPARAMSLICE  Slice scalars or row vectors for the provided indices.

    if isscalar(param)
        slice = repmat(param, size(idx));
    else
        slice = param(idx);
    end
end
