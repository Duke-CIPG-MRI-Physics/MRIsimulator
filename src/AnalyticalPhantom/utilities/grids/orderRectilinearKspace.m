function [kOrdered, tSamp] = orderRectilinearKspace(matrix_size, freq_phase_slice, dt, TR)
% orderRectilinearKspace  Generate rectilinear k-space sampling order.
%
%   [kOrdered, tSamp] = orderRectilinearKspace(matrix_size, freq_phase_slice, dt, TR)
%
%   Inputs:
%       matrix_size       : 1x3 vector [Nx, Ny, Nz] describing samples per
%                           k-space dimension (integer, > 0)
%       freq_phase_slice  : 1x3 permutation of [1 2 3] indicating which
%                           dimension is frequency-encode, phase-encode,
%                           and slice-encode, respectively
%       dt                : dwell time between consecutive frequency
%                           samples (seconds, > 0)
%       TR                : repetition time between the starts of
%                           successive frequency-encoding lines
%                           (seconds, > 0)
%
%   Outputs:
%       kOrdered : 3xN matrix of kx, ky, kz sample locations (1-based index units)
%       tSamp    : 1xN vector of cumulative sampling times for each sample
%
%   Notes:
%       - The function loops through the slice and phase dimensions in the
%         specified order, sampling the frequency dimension fastest. The
%         first sample occurs at t = 0.
%       - Each frequency sample within a line is separated by dt. The start
%         time of each subsequent frequency-encoding line (new phase or
%         slice position) is offset by TR.
%       - Output indices are simple 1..N indices along each dimension rather
%         than centered about zero.
%
%   Example:
%       [kOrdered, tSamp] = orderRectilinearKspace([4 3 2], [1 2 3], 1e-3, 5e-3);
%
    arguments
        matrix_size      (1,3) double {mustBeInteger, mustBePositive}
        freq_phase_slice (1,3) double {mustBeInteger, mustBePositive}
        dt               (1,1) double {mustBePositive}
        TR               (1,1) double {mustBePositive}
    end

    % Validate encoding order
    if ~isequal(sort(freq_phase_slice), [1 2 3])
        error('freq_phase_slice must be a permutation of [1 2 3].');
    end

    % Unpack encoding roles
    freqDim  = freq_phase_slice(1);
    phaseDim = freq_phase_slice(2);
    sliceDim = freq_phase_slice(3);

    % Sample counts along each encoded dimension
    nFreq  = matrix_size(freqDim);
    nPhase = matrix_size(phaseDim);
    nSlice = matrix_size(sliceDim);

    % Build an ordering with frequency varying fastest, then phase, then slice
    [phaseGrid, freqGrid, sliceGrid] = meshgrid(1:nPhase, 1:nFreq, 1:nSlice);
    encodedOrder = [freqGrid(:)'; phaseGrid(:)'; sliceGrid(:)'];

    % Map encoded order back to kx/ky/kz positions
    numSamples = numel(freqGrid);
    kOrdered   = zeros(3, numSamples);
    for dim = 1:3
        kOrdered(freq_phase_slice(dim), :) = encodedOrder(dim, :);
    end

    % Compute sampling timestamps matching the ordering (frequency samples
    % progress with dt; lines start every TR)
    nLines = numSamples / nFreq;
    lineStartTimes = repelem(0:nLines-1, nFreq) * TR;
    freqOffsets    = repmat(0:nFreq-1, 1, nLines) * dt;
    tSamp          = lineStartTimes + freqOffsets;
end
