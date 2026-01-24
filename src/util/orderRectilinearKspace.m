function [tSamp] = orderRectilinearKspace(matrix_size_fps, dt, TR)
% orderRectilinearKspace  Generate rectilinear k-space sampling order.
%
%   [kOrdered, tSamp] = orderRectilinearKspace(matrix_size, freq_phase_slice, dt, TR)
%
%   Inputs:
%       matrix_size_fps   : 1x3 vector [Nfreq, Nphase, Nslice] describing samples per
%                           k-space dimension (integer, > 0)
%       dt                : dwell time between consecutive frequency
%                           samples (seconds, > 0)
%       TR                : repetition time between the starts of
%                           successive frequency-encoding lines
%                           (seconds, > 0)
%
%   Outputs:
%       encodedOrder_fps : 3xN matrix of kfreq, kphase, kslice sample locations
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
        matrix_size_fps      (1,3) double {mustBeInteger, mustBePositive}
        dt               (1,1) double {mustBePositive}
        TR               (1,1) double {mustBePositive}
    end

    % Compute sampling timestamps matching the ordering (frequency samples
    % progress with dt; lines start every TR)
    nFreq = matrix_size_fps(1);
    numSamples = prod(matrix_size_fps);
    nLines = numSamples / nFreq;
    lineStartTimes = repelem(0:nLines-1, nFreq) * TR;
    freqOffsets    = repmat(0:nFreq-1, 1, nLines) * dt;
    tSamp          = lineStartTimes + freqOffsets;
end
