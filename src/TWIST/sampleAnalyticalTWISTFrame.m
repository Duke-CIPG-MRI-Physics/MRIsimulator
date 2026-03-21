function frameSamples = sampleAnalyticalTWISTFrame( ...
    phantom, Sampling_Table, twistPlan, frameIndex, ...
    kFreq, kPhase, kSlice, freqPhaseSlice, noiseSigma, maxChunkSize)
% sampleAnalyticalTWISTFrame  Evaluate one TWIST frame on demand.
%   frameSamples = sampleAnalyticalTWISTFrame(phantom, Sampling_Table, twistPlan, frameIndex,
%       kFreq, kPhase, kSlice, freqPhaseSlice, noiseSigma, maxChunkSize)
%   computes the raw complex samples for a single frame and splits them
%   into region A and the present B subsets.

arguments
    phantom
    Sampling_Table table
    twistPlan (1,1) struct
    frameIndex (1,1) {mustBePositive, mustBeInteger}
    kFreq
    kPhase
    kSlice
    freqPhaseSlice (1,3) {mustBeNumeric, mustBeInteger, mustBePositive}
    noiseSigma (1,1) {mustBeNumeric, mustBeNonnegative}
    maxChunkSize (1,1) {mustBeNumeric, mustBePositive, mustBeInteger}
end

frameRows = twistPlan.frameRowStart(frameIndex):twistPlan.frameRowStop(frameIndex);
nSamples = numel(frameRows);

kXYZ = zeros(3, nSamples, 'like', kFreq);
kXYZ(freqPhaseSlice(1), :) = kFreq(Sampling_Table.Frequency(frameRows));
kXYZ(freqPhaseSlice(2), :) = kPhase(Sampling_Table.("Row (phase)")(frameRows));
kXYZ(freqPhaseSlice(3), :) = kSlice(Sampling_Table.("Column (slice)")(frameRows));

t_s = Sampling_Table.Timing(frameRows)';
if isa(kFreq, 'gpuArray')
    t_s = gpuArray(t_s);
end

sampledValues = phantom.kspaceAtTime( ...
    kXYZ(1, :), ...
    kXYZ(2, :), ...
    kXYZ(3, :), ...
    t_s, ...
    maxChunkSize)';

if noiseSigma > 0
    noise = noiseSigma * ( ...
        randn(size(sampledValues), 'like', sampledValues) + ...
        1i * randn(size(sampledValues), 'like', sampledValues));
    sampledValues = sampledValues + noise;
end

localBj = Sampling_Table.Bj(frameRows);
frameSamples = struct("aValues", sampledValues(find(localBj == 0)), ...
    "bValues", {cell(1, twistPlan.nBSubsets)});

presentSubsets = unique(localBj(localBj > 0))';
for iSubset = presentSubsets
    frameSamples.bValues{iSubset} = sampledValues(find(localBj == iSubset));
end
end
