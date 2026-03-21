function twistPlan = prepareTWISTViewSharingPlan(Sampling_Table, Matrix_Size_Acquired, shareMode, tieBreaker)
% prepareTWISTViewSharingPlan  Build lightweight TWIST view-sharing metadata.
%   twistPlan = prepareTWISTViewSharingPlan(Sampling_Table, Matrix_Size_Acquired, shareMode)
%   returns compact per-frame metadata for streaming TWIST reconstruction.
%
%   Inputs
%   ------
%   Sampling_Table : table
%       TWIST sampling table after frequency encoding and timing assembly.
%   Matrix_Size_Acquired : [1x3] integer
%       Acquired matrix size in [frequency, phase, slice].
%   shareMode : string
%       "forward", "reverse", or "symmetric".
%   tieBreaker : string
%       Used only for symmetric mode when two candidate source frames are
%       equally close in time. Supported values: "future", "past".
%
%   Output
%   ------
%   twistPlan.frameNumbers          : [1xN] frame labels from Sampling_Table
%   twistPlan.frameRowStart         : [1xN] starting table row for each frame
%   twistPlan.frameRowStop          : [1xN] ending table row for each frame
%   twistPlan.frameCenterTimes_s    : [1xN] center-of-k-space sample time
%   twistPlan.frameIsFull           : [1xN] true when all B subsets are present
%   twistPlan.acquiredSubsetByFrame : [1xN] current Bj for partial frames, 0 otherwise
%   twistPlan.aLinearIdx            : linear indices for region A
%   twistPlan.bLinearIdxBySubset    : linear indices for each B subset
%   twistPlan.sourceFrameMap        : source frame for each frame/subset pair
%   twistPlan.readyFrameByFrame     : earliest acquired frame after which a
%                                     symmetric frame can be reconstructed

arguments
    Sampling_Table table
    Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}
    shareMode (1,1) string
    tieBreaker (1,1) string = "future"
end

shareMode = lower(shareMode);
tieBreaker = lower(tieBreaker);

validShareModes = ["forward", "reverse", "symmetric"];
if ~ismember(shareMode, validShareModes)
    error("TWIST:prepareTWISTViewSharingPlan:InvalidShareMode", ...
        "shareMode must be one of %s.", char(strjoin(validShareModes, ", ")));
end

validTieBreakers = ["future", "past"];
if ~ismember(tieBreaker, validTieBreakers)
    error("TWIST:prepareTWISTViewSharingPlan:InvalidTieBreaker", ...
        "tieBreaker must be one of %s.", char(strjoin(validTieBreakers, ", ")));
end

frameNumbers = unique(Sampling_Table.Frame, "stable")';
nFrames = numel(frameNumbers);
nBSubsets = max(Sampling_Table.Bj);

frameChangeRows = [1; find(diff(Sampling_Table.Frame) ~= 0) + 1];
frameRowStart = frameChangeRows';
frameRowStop = [frameChangeRows(2:end) - 1; height(Sampling_Table)]';

frameCenterTimes_s = nan(1, nFrames);
frameIsFull = false(1, nFrames);
acquiredSubsetByFrame = zeros(1, nFrames);
frameHasSubset = false(nFrames, nBSubsets);

kspaceCenter = floor(Matrix_Size_Acquired / 2) + 1;
aLinearIdx = [];
bLinearIdxBySubset = cell(1, nBSubsets);

for iFrame = 1:nFrames
    frameRows = frameRowStart(iFrame):frameRowStop(iFrame);
    localBj = Sampling_Table.Bj(frameRows);

    if isempty(aLinearIdx)
        aMask = (localBj == 0);
        aLinearIdx = sub2ind(Matrix_Size_Acquired, ...
            Sampling_Table.Frequency(frameRows(aMask)), ...
            Sampling_Table.("Row (phase)")(frameRows(aMask)), ...
            Sampling_Table.("Column (slice)")(frameRows(aMask)));
    end

    presentSubsets = unique(localBj(localBj > 0))';
    frameIsFull(iFrame) = (nBSubsets == 0) || (numel(presentSubsets) == nBSubsets);
    if numel(presentSubsets) == 1
        acquiredSubsetByFrame(iFrame) = presentSubsets;
    end

    for iSubset = 1:nBSubsets
        subsetMask = (localBj == iSubset);
        frameHasSubset(iFrame, iSubset) = any(subsetMask);

        if isempty(bLinearIdxBySubset{iSubset}) && any(subsetMask)
            bLinearIdxBySubset{iSubset} = sub2ind(Matrix_Size_Acquired, ...
                Sampling_Table.Frequency(frameRows(subsetMask)), ...
                Sampling_Table.("Row (phase)")(frameRows(subsetMask)), ...
                Sampling_Table.("Column (slice)")(frameRows(subsetMask)));
        end
    end

    centerMask = ...
        Sampling_Table.Frequency(frameRows) == kspaceCenter(1) & ...
        Sampling_Table.("Row (phase)")(frameRows) == kspaceCenter(2) & ...
        Sampling_Table.("Column (slice)")(frameRows) == kspaceCenter(3);
    if ~any(centerMask)
        error("TWIST:prepareTWISTViewSharingPlan:MissingCenterSample", ...
            "Frame %d does not contain the center of k-space.", frameNumbers(iFrame));
    end
    frameCenterTimes_s(iFrame) = Sampling_Table.Timing(frameRows(find(centerMask, 1, "first")));
end

sourceFrameMap = zeros(nFrames, nBSubsets);
for iFrame = 1:nFrames
    for iSubset = 1:nBSubsets
        candidateFrames = find(frameHasSubset(:, iSubset));
        if isempty(candidateFrames)
            error("TWIST:prepareTWISTViewSharingPlan:MissingSubsetSamples", ...
                "No samples were found for B%d.", iSubset);
        end

        switch shareMode
            case "forward"
                eligibleFrames = candidateFrames(candidateFrames <= iFrame);
                if isempty(eligibleFrames)
                    eligibleFrames = candidateFrames(1);
                end
                sourceFrameMap(iFrame, iSubset) = eligibleFrames(end);

            case "reverse"
                eligibleFrames = candidateFrames(candidateFrames >= iFrame);
                if isempty(eligibleFrames)
                    eligibleFrames = candidateFrames(end);
                end
                sourceFrameMap(iFrame, iSubset) = eligibleFrames(1);

            case "symmetric"
                candidateTimes = frameCenterTimes_s(candidateFrames);
                timeDelta = abs(candidateTimes - frameCenterTimes_s(iFrame));
                minDelta = min(timeDelta);
                tiedFrames = candidateFrames(timeDelta == minDelta);

                if numel(tiedFrames) == 1
                    sourceFrameMap(iFrame, iSubset) = tiedFrames;
                elseif tieBreaker == "future"
                    futureFrames = tiedFrames(tiedFrames >= iFrame);
                    if isempty(futureFrames)
                        sourceFrameMap(iFrame, iSubset) = tiedFrames(end);
                    else
                        sourceFrameMap(iFrame, iSubset) = futureFrames(1);
                    end
                else
                    pastFrames = tiedFrames(tiedFrames <= iFrame);
                    if isempty(pastFrames)
                        sourceFrameMap(iFrame, iSubset) = tiedFrames(1);
                    else
                        sourceFrameMap(iFrame, iSubset) = pastFrames(end);
                    end
                end
        end
    end
end

if nBSubsets == 0
    readyFrameByFrame = (1:nFrames)';
else
    readyFrameByFrame = max([(1:nFrames)', max(sourceFrameMap, [], 2)], [], 2);
end

twistPlan = struct;
twistPlan.frameNumbers = frameNumbers;
twistPlan.nFrames = nFrames;
twistPlan.nBSubsets = nBSubsets;
twistPlan.frameRowStart = frameRowStart;
twistPlan.frameRowStop = frameRowStop;
twistPlan.frameCenterTimes_s = frameCenterTimes_s;
twistPlan.frameIsFull = frameIsFull;
twistPlan.acquiredSubsetByFrame = acquiredSubsetByFrame;
twistPlan.aLinearIdx = aLinearIdx;
twistPlan.bLinearIdxBySubset = bLinearIdxBySubset;
twistPlan.sourceFrameMap = sourceFrameMap;
twistPlan.readyFrameByFrame = readyFrameByFrame;
end
