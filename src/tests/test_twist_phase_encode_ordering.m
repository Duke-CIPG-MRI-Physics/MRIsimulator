function test_twist_phase_encode_ordering()
% test_twist_phase_encode_ordering  Validate physical TWIST PE ordering.
%   test_twist_phase_encode_ordering() checks that the phase/slice ordering
%   is built on the physical k-space grid, that Partial Fourier mask counts
%   are correct, and that contiguous Bj assignment produces the expected
%   per-frame B trajectories.

    Matrix_Size_Acquired = [1, 12, 14];
    FOV_acquired = [1, 240, 140];
    pA = 0.20;
    pB = 0.25;
    R = [1, 1];
    PF_Factor = [1, 1];
    orderingOptions = getTWISTOrderingOptions(struct( ...
        'radialBinWidthMode', "max", ...
        'bSubsetAssignment', "contiguous"));

    [phaseEncodeTable, orderingInfo] = buildTWISTPhaseEncodeTable( ...
        Matrix_Size_Acquired, FOV_acquired, PF_Factor, R, orderingOptions);
    [regionA, ~] = getRegionA( ...
        Matrix_Size_Acquired, FOV_acquired, pA, PF_Factor, R, orderingOptions);
    twistTable = TWIST( ...
        pA, pB, Matrix_Size_Acquired, FOV_acquired, R, PF_Factor, orderingOptions);

    centerIndex = sub2ind(Matrix_Size_Acquired(2:3), ...
        floor(Matrix_Size_Acquired(2) / 2) + 1, ...
        floor(Matrix_Size_Acquired(3) / 2) + 1);
    assert(phaseEncodeTable.LinearIndex(1) == centerIndex, ...
        'The first ordered phase/slice point should be the DC sample.');
    assert(abs(phaseEncodeTable.SpatialFrequencyRadius(1)) < 1e-12, ...
        'The first ordered phase/slice point should have zero spatial frequency radius.');
    assert(abs(orderingInfo.radialBinWidth - max(1 / FOV_acquired(2), 1 / FOV_acquired(3))) < 1e-12, ...
        'The max shell-width mode should use the larger phase/slice k-space pixel size.');

    expectedAPoints = max(1, round(pA * Matrix_Size_Acquired(2) * Matrix_Size_Acquired(3)));
    assert(nnz(regionA) == expectedAPoints, ...
        'Region A should contain the requested number of central PE samples when no acceleration is used.');

    acceleratedTable = buildTWISTPhaseEncodeTable( ...
        [1, 20, 20], [1, 200, 200], [0.75, 0.75], [1, 1], orderingOptions);
    assert(nnz(acceleratedTable.IsKeptAfterAcceleration) == 225, ...
        'Partial Fourier should retain the expected number of phase/slice samples.');

    sortedBRows = phaseEncodeTable(~regionA(phaseEncodeTable.LinearIndex), :);
    nBSubsets = round(1 / pB);
    subsetEdges = round(linspace(0, height(sortedBRows), nBSubsets + 1));
    for iSubset = 1:nBSubsets
        currentChunk = sortedBRows(subsetEdges(iSubset) + 1:subsetEdges(iSubset + 1), :);
        expectedOrder = buildReferenceTrajectory(currentChunk, "inside_out");

        actualFrameRows = twistTable(twistTable.Frame == iSubset & twistTable.Bj > 0, :);
        assert(isequal(actualFrameRows.("Linear Index"), expectedOrder.LinearIndex), ...
            'Frame %d B ordering did not match the contiguous chunk trajectory.', iSubset);
    end
end

function orderedRows = buildReferenceTrajectory(sortedRows, traversalMode)
% buildReferenceTrajectory  Mirror the TWIST in/out trajectory logic.

    switch traversalMode
        case "outside_in"
            baseRows = flipud(sortedRows);
        case "inside_out"
            baseRows = sortedRows;
        otherwise
            error('Unsupported traversal mode.');
    end

    orderedRows = [baseRows(1:2:end, :); flipud(baseRows(2:2:end, :))];
end
