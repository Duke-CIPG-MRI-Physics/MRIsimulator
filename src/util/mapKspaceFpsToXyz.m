function [kXYZ, fpsToXyz] = mapKspaceFpsToXyz(kFPS, freqPhaseSlice)
% mapKspaceFpsToXyz Map k-space from FPS ordering to physical XYZ ordering.
%
%   [kXYZ, fpsToXyz] = mapKspaceFpsToXyz(kFPS, freqPhaseSlice)
%
%   Inputs:
%       kFPS          : [3xN] double, k-space samples in frequency, phase,
%                        and slice order (cycles/mm).
%       freqPhaseSlice: [1x3] integer, encoding directions for
%                        [frequency, phase, slice].
%                        1 = R/L, 2 = A/P, 3 = S/I.
%
%   Outputs:
%       kXYZ     : [3xN] double, k-space samples in physical XYZ ordering.
%       fpsToXyz : [1x3] integer, permutation mapping FPS -> XYZ.
%
%   Example:
%       [kXYZ, fpsToXyz] = mapKspaceFpsToXyz(kFPS, [2 1 3]);

    arguments
        kFPS (3,:) double
        freqPhaseSlice (1,3) double {mustBeInteger, mustBePositive}
    end

    validateattributes(freqPhaseSlice, {'double'}, ...
        {'integer', '>=', 1, '<=', 3, 'numel', 3}, mfilename, 'freqPhaseSlice');

    kXYZ = zeros(3, size(kFPS, 2));
    kXYZ(freqPhaseSlice(1), :) = kFPS(1, :);
    kXYZ(freqPhaseSlice(2), :) = kFPS(2, :);
    kXYZ(freqPhaseSlice(3), :) = kFPS(3, :);

    fpsToXyz = zeros(1, 3);
    fpsToXyz(freqPhaseSlice) = 1:3;
end
