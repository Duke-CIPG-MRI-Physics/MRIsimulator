function [regionA, frequency_table] = getRegionA( ...
    Matrix_Size_Acquired, FOV_acquired, pA, PF_Factor, R, radialBinWidthMode)
% getRegionA  Calculate the TWIST region-A mask on the phase/slice plane.
%   [regionA, frequency_table] = getRegionA(Matrix_Size_Acquired, ...
%       FOV_acquired, pA, PF_Factor, R, radialBinWidthMode) returns a logical
%   mask for region A and the underlying phase/slice ordering table.
%
%   Region A is selected from the lowest-radius entries on the physical
%   phase/slice k-space grid after quantizing radii into integer shells.
%   The number of A points is still controlled by pA, measured relative to
%   the number of samples expected after GRAPPA and Partial Fourier masks.

    arguments
        Matrix_Size_Acquired (1,3) {mustBeNumeric, mustBePositive, mustBeInteger}
        FOV_acquired (1,3) {mustBeNumeric, mustBePositive}
        pA (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(pA, 0), mustBeLessThanOrEqual(pA, 1)}
        PF_Factor (1,2) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(PF_Factor, 1), mustBeGreaterThan(PF_Factor, 0.5)}
        R (1,2) {mustBeNumeric, mustBePositive, mustBeInteger}
        radialBinWidthMode (1,1) string = "max"
    end

    [frequency_table, ~] = buildTWISTPhaseEncodeTable( ...
        Matrix_Size_Acquired, FOV_acquired, PF_Factor, R, radialBinWidthMode);

    retainedRows = find(frequency_table.IsKeptAfterAcceleration);
    n_pixels_acquired = numel(retainedRows);
    n_pixels_in_A = min(n_pixels_acquired, max(1, round(pA * n_pixels_acquired)));

    regionA = false(Matrix_Size_Acquired(2), Matrix_Size_Acquired(3));
    regionA(frequency_table.LinearIndex(retainedRows(1:n_pixels_in_A))) = true;
end
