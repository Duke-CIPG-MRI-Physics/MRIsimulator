function test_cardiac_waveform_chunk_consistency()
%TEST_CARDIAC_WAVEFORM_CHUNK_CONSISTENCY  Ensure chunked cardiac calls match full call.
%
%   test_cardiac_waveform_chunk_consistency()
%
%   Verifies equivalence between full time-vector evaluation and chunked
%   evaluation using absolute-time anchoring, and checks streaming phase
%   carryover using startPhase_cycles for constant-HR operation.

    %% Deterministic setup
    baseOpts = struct();
    baseOpts.HR_bpm = 75;
    baseOpts.EDV_ml = 150;
    baseOpts.ESV_ml = 75;
    baseOpts.systFrac = 0.35;
    baseOpts.q_ED = 50/27;
    baseOpts.GLS_peak = -0.20;
    baseOpts.GCS_peak = -0.25;

    t_s = 0.2:0.05:5.2;

    %% Full run reference
    fullParams = cardiac_ellipsoid_waveform(t_s, baseOpts);

    %% Chunked with absolute time (default reference)
    splitIndices = [1, 13; 14, 37; 38, numel(t_s)];
    chunkedAbsolute = assembleChunkedResult(t_s, splitIndices, baseOpts);

    assert(max(abs(fullParams.a_mm(:) - chunkedAbsolute.a_mm(:))) < 1e-12, ...
        'Absolute-time chunking did not match full run for a_mm.');
    assert(max(abs(fullParams.b_mm(:) - chunkedAbsolute.b_mm(:))) < 1e-12, ...
        'Absolute-time chunking did not match full run for b_mm.');
    assert(max(abs(fullParams.c_mm(:) - chunkedAbsolute.c_mm(:))) < 1e-12, ...
        'Absolute-time chunking did not match full run for c_mm.');

    %% Chunked with explicit state carryover (startPhase_cycles)
    chunkedState = assembleChunkedResultWithState(t_s, splitIndices, baseOpts);

    assert(max(abs(fullParams.a_mm(:) - chunkedState.a_mm(:))) < 1e-12, ...
        'State-carry chunking did not match full run for a_mm.');
    assert(max(abs(fullParams.b_mm(:) - chunkedState.b_mm(:))) < 1e-12, ...
        'State-carry chunking did not match full run for b_mm.');
    assert(max(abs(fullParams.c_mm(:) - chunkedState.c_mm(:))) < 1e-12, ...
        'State-carry chunking did not match full run for c_mm.');

    %% Single-sample evaluation no longer errors
    t_single_s = 0;
    singleParams = cardiac_ellipsoid_waveform(t_single_s, baseOpts);
    assert(isscalar(singleParams.a_mm) && isscalar(singleParams.b_mm) && isscalar(singleParams.c_mm), ...
        'Single-sample evaluation should return scalar shape parameters.');
end

function chunkedParams = assembleChunkedResult(t_s, splitIndices, opts)
%ASSEMBLECHUNKEDRESULT  Evaluate waveform in chunks and concatenate outputs.

    numSamples = numel(t_s);
    a_mm = zeros(1, numSamples);
    b_mm = zeros(1, numSamples);
    c_mm = zeros(1, numSamples);

    for idxChunk = 1:size(splitIndices, 1)
        idx = splitIndices(idxChunk, 1):splitIndices(idxChunk, 2);
        chunkParams = cardiac_ellipsoid_waveform(t_s(idx), opts);
        a_mm(idx) = chunkParams.a_mm;
        b_mm(idx) = chunkParams.b_mm;
        c_mm(idx) = chunkParams.c_mm;
    end

    chunkedParams = struct('a_mm', a_mm, 'b_mm', b_mm, 'c_mm', c_mm);
end

function chunkedParams = assembleChunkedResultWithState(t_s, splitIndices, opts)
%ASSEMBLECHUNKEDRESULTWITHSTATE  Evaluate waveform in chunks with phase carryover.

    numSamples = numel(t_s);
    a_mm = zeros(1, numSamples);
    b_mm = zeros(1, numSamples);
    c_mm = zeros(1, numSamples);

    frequency_Hz = opts.HR_bpm / 60;
    currentStartPhase_cycles = mod(frequency_Hz * t_s(splitIndices(1, 1)), 1);

    for idxChunk = 1:size(splitIndices, 1)
        idx = splitIndices(idxChunk, 1):splitIndices(idxChunk, 2);

        chunkOpts = opts;
        chunkOpts.startPhase_cycles = currentStartPhase_cycles;

        chunkParams = cardiac_ellipsoid_waveform(t_s(idx), chunkOpts);
        a_mm(idx) = chunkParams.a_mm;
        b_mm(idx) = chunkParams.b_mm;
        c_mm(idx) = chunkParams.c_mm;

        if idxChunk < size(splitIndices, 1)
            nextStart = splitIndices(idxChunk + 1, 1);
            deltaToNextStart_s = t_s(nextStart) - t_s(idx(1));
            currentStartPhase_cycles = mod(currentStartPhase_cycles + ...
                frequency_Hz * deltaToNextStart_s, 1);
        end
    end

    chunkedParams = struct('a_mm', a_mm, 'b_mm', b_mm, 'c_mm', c_mm);
end
