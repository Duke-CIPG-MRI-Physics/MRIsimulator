function [a_mm, b_mm, phase] = cardiac_ellipsoid_waveform(t_s, opts)
%CARDIAC_ELLIPSOID_WAVEFORM  LV ellipsoid semi-axes driven by volume + strain.
%
%   [a_mm, b_mm, phase] = cardiac_ellipsoid_waveform(t_s, opts)
%
%   Inputs:
%       t_s   - time [s], strictly increasing (length N)
%       opts  - struct with required fields HR_bpm, EDV_ml, ESV_ml. Each can be
%               a scalar, vector (length N), or matrix with N columns. For
%               matrices, rows correspond to independent cardiac waveforms. The
%               following optional fields set model parameters:
%                   .systFrac  - systolic fraction of R-R interval (0<beta<1),
%                                default 0.35
%                   .q_ED      - ED aspect ratio a0/b0 (long/short), default 2.5
%                   .GLS_peak  - peak global longitudinal strain (negative),
%                                default -0.20
%                   .GCS_peak  - peak global circumferential strain (negative),
%                                default -0.25
%
%   Outputs (matching the number of waveform rows in opts inputs):
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
        t_s (1,:) double {mustBeReal, mustBeFinite}
        opts struct = struct()
    end

    requiredFields = {'HR_bpm', 'EDV_ml', 'ESV_ml'};
    for idxField = 1:numel(requiredFields)
        if ~isfield(opts, requiredFields{idxField}) || isempty(opts.(requiredFields{idxField}))
            error('cardiac_ellipsoid_waveform:MissingField', ...
                'opts.%s is required and cannot be empty.', requiredFields{idxField});
        end
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

    HR_bpm = opts.HR_bpm;
    EDV_ml = opts.EDV_ml;
    ESV_ml = opts.ESV_ml;

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

    [HR_bpm, HR_rows] = validateLengthOrScalar(HR_bpm, numSamples, 'opts.HR_bpm', numSamples - 1);
    [EDV_ml, EDV_rows] = validateLengthOrScalar(EDV_ml, numSamples, 'opts.EDV_ml');
    [ESV_ml, ESV_rows] = validateLengthOrScalar(ESV_ml, numSamples, 'opts.ESV_ml');

    numChannels = max([1, HR_rows, EDV_rows, ESV_rows]);

    ensureChannelCompatibility(HR_rows, numChannels, 'opts.HR_bpm');
    ensureChannelCompatibility(EDV_rows, numChannels, 'opts.EDV_ml');
    ensureChannelCompatibility(ESV_rows, numChannels, 'opts.ESV_ml');

    EDV_full = extractParamSlice(EDV_ml, 1:numSamples, numChannels);
    ESV_full = extractParamSlice(ESV_ml, 1:numSamples, numChannels);
    if any(ESV_full >= EDV_full, 'all')
        error('ESV_ml must be < EDV_ml at all time points.');
    end

    % ---------------------------------------------------------------------
    % 1) Build cumulative phase from HR(t)
    % ---------------------------------------------------------------------
    f_Hz = HR_bpm / 60;     % [Hz]

    f_interval = buildIntervalFrequencies(f_Hz, numel(dt), numSamples, numChannels);

    phaseIncrement = 2*pi .* f_interval .* dt;
    phase = mod([zeros(numChannels, 1), cumsum(phaseIncrement, 2)], 2*pi) / (2*pi);
    clear dt f_interval phaseIncrement;
    
    % ---------------------------------------------------------------------
    % 2) Volume waveform via piecewise half-cosines
    % ---------------------------------------------------------------------
    % Chunk calculations to ~0.001 Gb of temporary data so memory impact is
    % minimal.
    chunkSize = 15625;
    V_ED = -Inf(numChannels, 1);

    for idxStart = 1:chunkSize:numSamples
        idxEnd = min(idxStart + chunkSize - 1, numSamples);
        idx = idxStart:idxEnd;

        waveformChunk = zeros(numChannels, numel(idx));
        eps_L_chunk = zeros(numChannels, numel(idx));
        eps_C_chunk = zeros(numChannels, numel(idx));

        [waveformChunk, eps_L_chunk, eps_C_chunk] = buildWaveformAndStrain( ...
            phase(:, idx), beta, GLS_peak, GCS_peak, waveformChunk, ...
            eps_L_chunk, eps_C_chunk);

        V_chunk = computeVolumeChunk(EDV_ml, ESV_ml, idx, waveformChunk, numChannels);

        V_chunk_max = max(V_chunk, [], 2);
        V_ED = max(V_ED, V_chunk_max);
    end

    % ---------------------------------------------------------------------
    % 3) Baseline ED ellipsoid from EDV and aspect ratio q = a0/b0
    % ---------------------------------------------------------------------

    % Prolate spheroid: a0 = q * b0, c0 = b0
    % V_ED = (4/3)*pi*a0*b0^2 = (4/3)*pi*q*b0^3  =>  solve for b0
    b0 = ((3*V_ED) / (4*pi*q)).^(1/3);
    a0 = q .* b0;
    % c0 = b0;  % implied

    % ---------------------------------------------------------------------
    % 4-5) Strain-driven axes + isotropic scale factor (chunked)
    % ---------------------------------------------------------------------
    a_mm = zeros(numChannels, numSamples);
    b_mm = zeros(numChannels, numSamples);

    for idxStart = 1:chunkSize:numSamples
        idxEnd = min(idxStart + chunkSize - 1, numSamples);
        idx = idxStart:idxEnd;

        waveformChunk = zeros(numChannels, numel(idx));
        eps_L_chunk = zeros(numChannels, numel(idx));
        eps_C_chunk = zeros(numChannels, numel(idx));

        [waveformChunk, eps_L_chunk, eps_C_chunk] = buildWaveformAndStrain( ...
            phase(:, idx), beta, GLS_peak, GCS_peak, waveformChunk, ...
            eps_L_chunk, eps_C_chunk);

        V_chunk = computeVolumeChunk(EDV_ml, ESV_ml, idx, waveformChunk, numChannels);

        a_hat_chunk = a0 .* (1 + eps_L_chunk);
        b_hat_chunk = b0 .* (1 + eps_C_chunk);
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
function [vec, numRows] = validateLengthOrScalar(vec, primaryLength, name, altLength)
%VALIDATELENGTHORSCALAR  Ensure input is scalar, vector, or 2D with time in columns.

    if nargin < 4
        altLength = [];
    end

    if isscalar(vec)
        numRows = 0;
        return;
    end

    validLengths = [primaryLength, altLength];
    validLengths = validLengths(validLengths > 0);

    if isvector(vec)
        vec = reshape(vec, 1, numel(vec));
        if ~ismember(numel(vec), validLengths)
            if isempty(altLength)
                error('%s must be scalar or have %d elements.', name, primaryLength);
            else
                error('%s must be scalar or have %d elements (or %d if specified per interval).', ...
                    name, primaryLength, primaryLength - 1);
            end
        end
        numRows = 1;
        return;
    end

    if ndims(vec) > 2
        error('%s must be scalar, a vector, or a 2D matrix with time along one dimension.', name);
    end

    sz = size(vec);
    if ismember(sz(2), validLengths)
        numRows = sz(1);
    elseif ismember(sz(1), validLengths) && sz(2) ~= sz(1)
        vec = vec.';
        numRows = size(vec, 1);
    else
        error('%s must have %d columns (or rows) matching t_s.', name, primaryLength);
    end
end

% -------------------------------------------------------------------------
function ensureChannelCompatibility(numRows, numChannels, name)
%ENSURECHANNELCOMPATIBILITY  Verify channels align or can be broadcast.

    if numRows > 1 && numRows ~= numChannels
        error('%s must have %d rows (channels) or be scalar/1-row for broadcasting.', ...
            name, numChannels);
    end
end

% -------------------------------------------------------------------------
function f_interval = buildIntervalFrequencies(f_Hz, numIntervals, numSamples, numChannels)
%BUILDINTERVALFREQUENCIES  Align heart rate waveform with interval durations.

    if isscalar(f_Hz)
        f_interval = f_Hz * ones(numChannels, numIntervals);
        return;
    end

    hrRows = size(f_Hz, 1);
    numSamplesProvided = size(f_Hz, 2);
    if numSamplesProvided == numIntervals
        f_interval = f_Hz;
    elseif numSamplesProvided == numSamples
        f_interval = f_Hz(:, 1:end-1);
    else
        error(['HR_bpm must be scalar or have the same length as t_s or diff(t_s). ', ...
            'Received %d samples.'], numSamplesProvided);
    end

    if hrRows == 1 && numChannels > 1
        f_interval = repmat(f_interval, numChannels, 1);
    elseif hrRows > 1 && hrRows ~= numChannels
        error('opts.HR_bpm must have %d rows (channels) or be scalar/1-row for broadcasting.', ...
            numChannels);
    end
end

% -------------------------------------------------------------------------
function V_chunk = computeVolumeChunk(EDV_ml, ESV_ml, idx, waveformChunk, numChannels)
%COMPUTEVOLUMECHUNK  Volume waveform in cubic meters for a chunk of indices.

    EDV_slice = extractParamSlice(EDV_ml, idx, numChannels);
    ESV_slice = extractParamSlice(ESV_ml, idx, numChannels);
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

    if any(systMask, 'all')
        f_syst = 0.5*(1 - cos(pi*phaseChunk(systMask) / beta));              % 0 -> 1
        waveformChunk(systMask) = 0.5*(1 + cos(pi*phaseChunk(systMask) / beta));
        eps_L_chunk(systMask) = GLS_peak * f_syst;
        eps_C_chunk(systMask) = GCS_peak * f_syst;
    end

    if any(diasMask, 'all')
        shiftedPhase = phaseChunk(diasMask) - beta;
        f_dias = 0.5*(1 + cos(pi * shiftedPhase / (1 - beta)));                   % 1 -> 0
        waveformChunk(diasMask) = 0.5*(1 - cos(pi * shiftedPhase / (1 - beta)));
        eps_L_chunk(diasMask) = GLS_peak * f_dias;
        eps_C_chunk(diasMask) = GCS_peak * f_dias;
    end
end

% -------------------------------------------------------------------------
function slice = extractParamSlice(param, idx, numChannels)
%EXTRACTPARAMSLICE  Slice scalars or matrices for the provided indices.

    if isscalar(param)
        slice = repmat(param, numChannels, numel(idx));
        return;
    end

    paramRows = size(param, 1);
    baseSlice = param(:, idx);

    if paramRows == 1 && numChannels > 1
        slice = repmat(baseSlice, numChannels, 1);
    elseif paramRows == numChannels
        slice = baseSlice;
    else
        error('Parameter must have %d rows (channels) or be scalar/1-row for broadcasting.', ...
            numChannels);
    end
end
