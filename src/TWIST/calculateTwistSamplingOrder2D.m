function [A1_idx_outerIn, A2_idx_innerOut] = calculateTwistSamplingOrder2D(pA, pB, FOV_acquired, matrix_size_acquired, piAcceleration)

% Compute total number of phase/slice encodes
nPhaseSliceEncodes = prod(matrix_size_acquired(2:3));  % nPhase * nSlice

% The A region of TWIST considers parallel imaging acceleration, but
% oddly does not undersample the A region
totalAcceleration = prod(piAcceleration);
pA_eff = pA / totalAcceleration;
nVoxA = ceil(pA_eff*nPhaseSliceEncodes);

% Compute the spatial frequencies in the phase and slice directions in
% units of mm^-1
spatialFreq_phase = computeKspaceGrid1D(FOV_acquired(2), matrix_size_acquired(2));
spatialFreq_slice = computeKspaceGrid1D(FOV_acquired(3), matrix_size_acquired(3));


% Calculate the smallest step size
minStepSize = min(spatialFreq_phase(2)-spatialFreq_phase(1),...
    spatialFreq_slice(2)-spatialFreq_slice(1));

% Calculate mesh of k-space pixels according to minimum size
[mSlice, mPhase] = meshgrid(spatialFreq_slice/minStepSize,spatialFreq_phase/minStepSize);

%Calculate rounded R, theta for each point
[theta, R] = cart2pol(mPhase,mSlice);
R = round(R); % bins into rings of thickness 1

% Sort into rings of R, then sort by angle theta
[tempSorted, sortedIdx] = sortrows([R(:), theta(:)]);
tempSorted

%% Calculate A region pixels, which loop around the rings, first going from the outer edge of k-space
% into the center of k-space, then going from the center back out to
% the edge in an interleaved fashion. Note that the A region does not undersample, though the number
% pixels in A is affected by parallel imaging acceleration
A1_idx_outerIn = sortedIdx(1:2:nVoxA);
A1_idx_outerIn = flipud(A1_idx_outerIn);
A2_idx_innerOut = sortedIdx(2:2:nVoxA);

A_mask = zeros(matrix_size_acquired(2:3));
A_mask(A1_idx_outerIn) = 1;
A_mask(A2_idx_innerOut) = 2;

%% Calculate parallel imaging sampling
pi_phases = (abs(mod(mPhase,piAcceleration(1))) < 1E-10);
pi_slices = (abs(mod(mSlice,piAcceleration(2))) < 1E-10);
pi_samples = pi_phases & pi_slices;

figure();
subplot(1,3,1);
imagesc(A_mask)
subplot(1,3,2);
imagesc(pi_samples)


%% Calculate B region
% Frames covering B given pB
nFrames = max(1, ceil(1 / pB));

% B region the parallel imaging samples that are not in A
B_mask = pi_samples & (~A_mask)



nVoxB_all = nPhaseSliceEncodes - nVoxA;

% Approximate B encodes per frame (integer)
nVoxB_perFrame = ceil(nVoxB_all / nFrames);  % Note that some TWIST frames will have fewer voxels



bRegion = true(matrix_size_acquired(2:3));
bRegion(A1_idx_outerIn) = false;
bRegion(A2_idx_innerOut) = false;
bRegion(~pi_samples) = false();

% Sort B region in rings
[~, sortedBIdx] = sortrows([R(bRegion), theta(bRegion)]);

bSortedIdx = sortedIdx(nVoxA+1:end)
end