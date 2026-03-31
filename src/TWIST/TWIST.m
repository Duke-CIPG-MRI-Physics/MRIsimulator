function [TWIST_sampling_order] = TWIST( ...
    pA, pB, Matrix_Size_Acquired, FOV_acquired, R, PF_Factor, orderingOptions)
% TWIST  Generate phase/slice TWIST sampling order on the physical k-space grid.
%   TWIST_sampling_order = TWIST(pA, pB, Matrix_Size_Acquired, FOV_acquired,
%       R, PF_Factor, orderingOptions) returns the phase/slice sampling table
%   for one complete TWIST cycle: an initial full acquisition followed by
%   interleaved A + Bj partial frames.
%
%   The phase/slice encodes are ordered using physical spatial frequencies
%   [cycles / distance], not raw pixel offsets. Radii are quantized into
%   integer shells before angular sorting so anisotropic FOV and matrix size
%   affect the ordering correctly.

    arguments
        pA (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(pA, 0), mustBeLessThanOrEqual(pA, 1)}
        pB {mustBeMember(pB, [0, .1, .25, .33, .5])}
        Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}
        FOV_acquired (1,3) {mustBeNumeric, mustBePositive}
        R (1,2) {mustBeNumeric, mustBeInteger, mustBePositive}
        PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor, 1), mustBeGreaterThan(PF_Factor, .5)}
        orderingOptions (1,1) struct = struct()
    end

    if pB == 0
        Nb = 0;
    else
        Nb = round(1 / pB);
    end

    orderingOptions = getTWISTOrderingOptions(orderingOptions);

    %% Build the physical phase/slice ordering table
    [regionA, frequency_table] = getRegionA( ...
        Matrix_Size_Acquired, FOV_acquired, pA, PF_Factor, R, orderingOptions);

    sortedData = frequency_table(:, ...
        {'LinearIndex', 'SpatialFrequencyRadius', 'RadialBin', 'Theta', 'PhaseRow', 'SliceColumn'});
    sortedData.Properties.VariableNames = { ...
        'Linear Index', ...
        'SpatialFrequencyRadius', ...
        'RadialBin', ...
        'Theta', ...
        'Row (phase)', ...
        'Column (slice)'};
    sortedData.("Region A?") = regionA(sortedData.("Linear Index"));
    sortedData.Bj = zeros(height(sortedData), 1);

    sortedData_regionA = sortedData(sortedData.("Region A?"), :);
    sortedData_regionB = sortedData(~sortedData.("Region A?"), :);

    if Nb > 0 && ~isempty(sortedData_regionB)
        sortedData_regionB.Bj = assignBSubsets( ...
            height(sortedData_regionB), Nb, orderingOptions.bSubsetAssignment);
        sortedData.Bj(~sortedData.("Region A?")) = sortedData_regionB.Bj;
    end

    %% Build A and B trajectories
    kspaceSamplingOrder_A = buildBidirectionalTrajectory(sortedData_regionA, "outside_in");

    kspaceSamplingOrder_all_frames = [ ...
        sortedData([], :), ...
        table(zeros(0, 1), 'VariableNames', {'Frame'})];
    kspaceSamplingOrder_B = sortedData([], :);

    if Nb == 0
        kspaceSamplingOrder_B = buildBidirectionalTrajectory(sortedData_regionB, "inside_out");
    else
        for ii = 1:Nb
            B_current_Bj = sortedData_regionB(sortedData_regionB.Bj == ii, :);
            kspaceSamplingOrder_B_current_Bj = buildBidirectionalTrajectory(B_current_Bj, "inside_out");

            kspaceSamplingOrder_B = [kspaceSamplingOrder_B; kspaceSamplingOrder_B_current_Bj]; %#ok<AGROW>

            kspaceSamplingOrder_current_frame = [kspaceSamplingOrder_A; kspaceSamplingOrder_B_current_Bj]; %#ok<AGROW>
            kspaceSamplingOrder_current_frame.Frame = ii * ones(height(kspaceSamplingOrder_current_frame), 1);
            kspaceSamplingOrder_all_frames = [kspaceSamplingOrder_all_frames; kspaceSamplingOrder_current_frame]; %#ok<AGROW>
        end
    end

    %% Initial full acquisition anchor
    kspaceSamplingOrder_initial = [kspaceSamplingOrder_B; kspaceSamplingOrder_A];
    kspaceSamplingOrder_initial.Frame = zeros(height(kspaceSamplingOrder_initial), 1);

    kspaceSamplingOrder_full = [kspaceSamplingOrder_initial; kspaceSamplingOrder_all_frames];

    %% Cleanup output table
    TWIST_sampling_order = kspaceSamplingOrder_full(:, {'Linear Index', 'Bj', 'Frame'});
    [row, col] = ind2sub(Matrix_Size_Acquired(2:3), TWIST_sampling_order.("Linear Index"));
    TWIST_sampling_order.("Row (phase)") = row;
    TWIST_sampling_order.("Column (slice)") = col;
end

function bjSequence = assignBSubsets(numRegionBPoints, nBSubsets, assignmentMode)
% assignBSubsets  Split the sorted B list into Bj subsets.

    if numRegionBPoints < nBSubsets
        error("TWIST:InsufficientBPoints", ...
            "Region B contains %d points, which is fewer than the %d Bj subsets.", ...
            numRegionBPoints, nBSubsets);
    end

    switch assignmentMode
        case "contiguous"
            bjSequence = zeros(numRegionBPoints, 1);
            subsetEdges = round(linspace(0, numRegionBPoints, nBSubsets + 1));
            for iSubset = 1:nBSubsets
                firstIdx = subsetEdges(iSubset) + 1;
                lastIdx = subsetEdges(iSubset + 1);
                if firstIdx <= lastIdx
                    bjSequence(firstIdx:lastIdx) = iSubset;
                end
            end

        case "interleaved"
            bjSequence = mod((0:numRegionBPoints - 1)', nBSubsets) + 1;
    end
end

function trajectoryTable = buildBidirectionalTrajectory(sortedRows, traversalMode)
% buildBidirectionalTrajectory  Build one in/out trajectory from a sorted list.

    if isempty(sortedRows)
        trajectoryTable = sortedRows;
        return
    end

    switch traversalMode
        case "outside_in"
            orderedRows = flipud(sortedRows);
        case "inside_out"
            orderedRows = sortedRows;
        otherwise
            error("TWIST:InvalidTraversalMode", ...
                "Unknown traversal mode '%s'.", traversalMode);
    end

    firstLeg = orderedRows(1:2:end, :);
    secondLeg = flipud(orderedRows(2:2:end, :));
    trajectoryTable = [firstLeg; secondLeg];
end
