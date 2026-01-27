function [nVoxA, nVoxB_all, nVoxB_perFrame, nFrames] = calculateTwistRegionPixelCounts( ...
    pA, pB, acquiredMatrixSize, totalAcceleration)
%CALCULATETWISTREGIONPIXELCOUNTS Estimate TWIST A/B region sample counts per frame.
%
%   [nVoxA, nVoxB_all, nVoxB_perFrame, nFrames] = CALCULATETWISTREGIONPIXELCOUNTS( ...
%       pA, pB, acquiredMatrixSize, totalAcceleration)
%
%   Implements a simplified TWIST sampling model that partitions phase-slice
%   k-space encodes into:
%     - Region A: a central fraction sampled every frame
%     - Region B: the remaining encodes, interleaved across frames
%
%   This function returns integer sample counts based on the acquired (non-zero-
%   padded) matrix size and the requested A/B fractions.
%
%   Inputs
%   ------
%   pA : scalar in (0,1]
%       Nominal fraction of phase-slice encodes assigned to region A.
%       NOTE: This implementation applies totalAcceleration to region A via:
%             pA_eff = pA / totalAcceleration
%       which matches the original code behavior. If your intent is that A is
%       *fully sampled* while only B is accelerated, you should NOT scale pA.
%
%   pB : scalar in (0,1]
%       Fraction of region B sampled per frame. Typically nFrames ~= 1/pB.
%       This implementation sets:
%           nFrames = ceil(1/pB)
%       so some frames may sample slightly less or more than pB of B.
%
%   acquiredMatrixSize : 1x3 (or 1x2) positive integers
%       Acquisition matrix size as [nRead nPhase nSlice] (preferred) or
%       [nRead nPhase]. Only phase and slice are used for counting encodes.
%
%   totalAcceleration : scalar positive
%       Effective acceleration in the phase/slice encoding domain from
%       parallel imaging and/or simultaneous multi-slice (SMS), etc.
%
%   Outputs
%   -------
%   nVoxA : integer >= 0
%       Number of phase-slice encodes in region A sampled every frame.
%
%   nVoxB_all : integer >= 0
%       Total number of phase-slice encodes in region B (across all frames).
%
%   nVoxB_perFrame : integer >= 0
%       Approximate number of B encodes sampled per frame (ceil division).
%       Note: Because of integer rounding, not every frame must have exactly
%       this count.
%
%   nFrames : integer >= 1
%       Number of frames needed to cover region B given pB.
%
%   Reference
%   ---------
%   Song T, Laine AF, Chen Q, et al. Optimal k-space sampling for dynamic
%   contrast-enhanced MRI with an application to MR renography.
%   Magn Reson Med. 2009;61(5):1242-1248. doi:10.1002/mrm.21901
%
%   Roberto Carrascosa, Duke University, June 2025

% -------------------- validation --------------------
arguments
    pA (1,1) double {mustBeFinite, mustBeGreaterThan(pA,0), mustBeLessThanOrEqual(pA,1)}
    pB (1,1) double {mustBeFinite, mustBeGreaterThan(pB,0), mustBeLessThanOrEqual(pB,1)}
    acquiredMatrixSize double {mustBeFinite, mustBeInteger, mustBePositive}
    totalAcceleration (1,1) double {mustBeFinite, mustBePositive}
end

if numel(acquiredMatrixSize) < 3
    error('acquiredMatrixSize must have at least 3 values [nRead nPhase nSlice].');
end

% Compute total number of phase/slice encodes
nPhaseSliceEncodes = prod(acquiredMatrixSize(2:3));  % nPhase * nSlice

% The A region of TWIST considers parallel imaging acceleration, but
% oddly does not undersample the A region
pA_eff = pA / totalAcceleration;
nVoxA = ceil(pA_eff*nPhaseSliceEncodes);

% B region is remainder
nVoxB_all = nPhaseSliceEncodes - nVoxA;

% Frames covering B given pB
nFrames = max(1, ceil(1 / pB));

% Approximate B encodes per frame (integer)
nVoxB_perFrame = ceil(nVoxB_all / nFrames);  % Note that some TWIST frames will have fewer voxels
end
