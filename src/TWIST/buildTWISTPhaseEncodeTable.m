function [phaseEncodeTable, orderingInfo] = buildTWISTPhaseEncodeTable( ...
    Matrix_Size_Acquired, FOV_acquired, PF_Factor, R, orderingOptions)
% buildTWISTPhaseEncodeTable  Build physical phase/slice k-space ordering data.
%   [phaseEncodeTable, orderingInfo] = buildTWISTPhaseEncodeTable( ...
%       Matrix_Size_Acquired, FOV_acquired, PF_Factor, R, orderingOptions)
%   returns the phase/slice encode plane as a table sorted by quantized
%   radial distance to DC and then by angular position.
%
%   Inputs
%   ------
%   Matrix_Size_Acquired : [1x3] integer
%       Acquired matrix size in [frequency, phase, slice].
%   FOV_acquired : [1x3] double
%       FOV in [frequency, phase, slice]. Units must be consistent.
%   PF_Factor : [1x2] double
%       Partial Fourier factors in [phase, slice].
%   R : [1x2] integer
%       GRAPPA acceleration in [phase, slice].
%   orderingOptions : struct
%       Normalized TWIST ordering options from getTWISTOrderingOptions().
%
%   Outputs
%   -------
%   phaseEncodeTable : table
%       Sorted table with columns:
%           LinearIndex
%           PhaseRow
%           SliceColumn
%           PhaseFrequency
%           SliceFrequency
%           SpatialFrequencyRadius
%           RadialBin
%           Theta
%           IsKeptAfterAcceleration
%   orderingInfo : struct
%       Shell width and axis metadata used to build the table.

    arguments
        Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}
        FOV_acquired (1,3) {mustBeNumeric, mustBePositive}
        PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor, 1), mustBeGreaterThan(PF_Factor, 0.5)}
        R (1,2) {mustBeNumeric, mustBePositive, mustBeInteger}
        orderingOptions (1,1) struct = struct()
    end

    orderingOptions = getTWISTOrderingOptions(orderingOptions);

    phaseAxis = computeKspaceGrid1D(FOV_acquired(2), Matrix_Size_Acquired(2));
    sliceAxis = computeKspaceGrid1D(FOV_acquired(3), Matrix_Size_Acquired(3));
    [phaseGrid, sliceGrid] = ndgrid(phaseAxis, sliceAxis);

    if numel(phaseAxis) == 1
        dkPhase = 1 / FOV_acquired(2);
    else
        dkPhase = abs(phaseAxis(2) - phaseAxis(1));
    end
    if numel(sliceAxis) == 1
        dkSlice = 1 / FOV_acquired(3);
    else
        dkSlice = abs(sliceAxis(2) - sliceAxis(1));
    end

    switch orderingOptions.radialBinWidthMode
        case "max"
            radialBinWidth = max(dkPhase, dkSlice);
        case "min"
            radialBinWidth = min(dkPhase, dkSlice);
    end

    spatialFrequencyRadius = hypot(phaseGrid, sliceGrid);
    radialBin = round(spatialFrequencyRadius ./ radialBinWidth);
    theta = mod(-atan2(phaseGrid, sliceGrid), 2 * pi);

    [phaseRows, sliceColumns] = ndgrid(1:Matrix_Size_Acquired(2), 1:Matrix_Size_Acquired(3));
    linearIndex = reshape(1:numel(spatialFrequencyRadius), size(spatialFrequencyRadius));
    isKeptAfterAcceleration = buildAccelerationMask(Matrix_Size_Acquired, PF_Factor, R);

    phaseEncodeTable = table( ...
        linearIndex(:), ...
        phaseRows(:), ...
        sliceColumns(:), ...
        phaseGrid(:), ...
        sliceGrid(:), ...
        spatialFrequencyRadius(:), ...
        radialBin(:), ...
        theta(:), ...
        isKeptAfterAcceleration(:), ...
        'VariableNames', { ...
            'LinearIndex', ...
            'PhaseRow', ...
            'SliceColumn', ...
            'PhaseFrequency', ...
            'SliceFrequency', ...
            'SpatialFrequencyRadius', ...
            'RadialBin', ...
            'Theta', ...
            'IsKeptAfterAcceleration'});

    phaseEncodeTable = sortrows(phaseEncodeTable, ...
        {'RadialBin', 'Theta', 'SpatialFrequencyRadius', 'PhaseRow', 'SliceColumn'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend', 'ascend'});

    orderingInfo = struct();
    orderingInfo.phaseAxis = phaseAxis;
    orderingInfo.sliceAxis = sliceAxis;
    orderingInfo.dkPhase = dkPhase;
    orderingInfo.dkSlice = dkSlice;
    orderingInfo.radialBinWidth = radialBinWidth;
    orderingInfo.radialBinWidthMode = orderingOptions.radialBinWidthMode;
    orderingInfo.bSubsetAssignment = orderingOptions.bSubsetAssignment;
end

function accelerationMask = buildAccelerationMask(Matrix_Size_Acquired, PF_Factor, R)
% buildAccelerationMask  Predict the rows retained after PF and GRAPPA masks.

    phaseCount = Matrix_Size_Acquired(2);
    sliceCount = Matrix_Size_Acquired(3);

    rowMask = true(phaseCount, 1);
    colMask = true(1, sliceCount);

    if any(R ~= 1)
        rowMask = false(phaseCount, 1);
        colMask = false(1, sliceCount);
        rowMask(R(1):R(1):phaseCount) = true;
        colMask(R(2):R(2):sliceCount) = true;
    end

    if any(PF_Factor ~= 1)
        pfRowsKept = round(phaseCount * PF_Factor(1));
        pfColsKept = round(sliceCount * PF_Factor(2));

        pfRowMask = false(phaseCount, 1);
        pfColMask = false(1, sliceCount);
        pfRowMask(1:pfRowsKept) = true;
        pfColMask((sliceCount - pfColsKept + 1):sliceCount) = true;

        rowMask = rowMask & pfRowMask;
        colMask = colMask & pfColMask;
    end

    accelerationMask = repmat(rowMask, 1, sliceCount) & repmat(colMask, phaseCount, 1);
end
