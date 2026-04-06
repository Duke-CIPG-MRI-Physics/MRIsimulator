function [TWIST_sampling_order] = TWIST( ...
    pA, pB, Matrix_Size_Acquired, FOV_acquired, R, PF_Factor, radialBinWidthMode)
%Roberto Carrascosa, Duke University, June 2025
% This function implements the TWIST sampling scheme for MRI
%
%This function outputs a table of points, in the order they should be
%sampled for TWIST. Note that TWIST operates only on 3D acquisitions.
%Therefore, actual use of this function will require further implemenetation
%of frequency encoding, based on the user's needs.
%
%The inputs of this function are:
% pA --- defines the size of the central region that is sampled for
%   every measurement
%
% pB --- where pB defines the fraction of the exteriorvof k-space which is
%   sampled every measurement. There is a limited set of allowed inputs
%
% Matrix_Size_Acquired --- vector of form: [#frequency, #phase (rows), #slice (columns)]
%
% FOV_acquired --- vector of form: [freq FOV, phase FOV, slice FOV]
%   units do not matter
%
% R --- GRAPPA acceleration factor of form [phase, slice]
%   can be skipped with value of 1
%
% PF_Factor --- Partial Fourier acceleration factor of form [phase fraction, slice fraction]
%   can be skipped with value of 1
%
% radialBinWidthMode --- optional string, either "max" or "min", used to
%   set the radial shell width in physical k-space units
%
%---------------------------------------
%
%The output of this function is a table of form:
% Linear Index, Bj, Frame, Row(phase), Column(slice)
%
% It is one complete TWIST sequence: full k-space, followed by interleaved
% A and B region in correct order.

    arguments
        pA (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(pA, 0), mustBeLessThanOrEqual(pA, 1)}
        pB {mustBeMember(pB, [0, .1, .25, .33, .5])}
        Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}
        FOV_acquired (1,3) {mustBeNumeric, mustBePositive}
        R (1,2) {mustBeNumeric, mustBeInteger, mustBePositive}
        PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor, 1), mustBeGreaterThan(PF_Factor, .5)}
        radialBinWidthMode (1,1) string = "max"
    end

    if pB == 0
        Nb = 0;
    else
        Nb = round(1 / pB);
    end

    %% Build the physical phase/slice ordering table
    [regionA, frequency_table] = getRegionA( ...
        Matrix_Size_Acquired, FOV_acquired, pA, PF_Factor, R, radialBinWidthMode);

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

    sortedData_regionA = sortedData(sortedData.("Region A?"),:);
    sortedData_regionB = sortedData(~sortedData.("Region A?"),:);

    if Nb > 0 && ~isempty(sortedData_regionB)
        sortedData_regionB.Bj = assignBSubsets(height(sortedData_regionB), Nb);
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

function bjSequence = assignBSubsets(numRegionBPoints, nBSubsets)
% assignBSubsets  Split the sorted B list into interleaved Bj subsets.

    if numRegionBPoints < nBSubsets
        error("TWIST:InsufficientBPoints", ...
            "Region B contains %d points, which is fewer than the %d Bj subsets.", ...
            numRegionBPoints, nBSubsets);
    end

    bjSequence = mod((0:numRegionBPoints - 1)', nBSubsets) + 1;
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
