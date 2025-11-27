function height_mm = computeContrastWashIn(t_s, radius_mm, V_contrast_mm3)
% computeContrastWashIn
%   Compute the vessel height (length) required to contain a time-varying
%   contrast volume, assuming a circular cylindrical cross-section.
%
%   height_mm = computeContrastWashIn(t_s, radius_mm, V_contrast_mm3)
%
%   Inputs:
%       t_s             : Time vector [s], numeric, real. Only the length is
%                         used for validation; values need not be monotonic.
%       radius_mm       : Vessel radius over time [mm], same length as t_s.
%       V_contrast_mm3  : Contrast volume over time [mm^3], same length as
%                         t_s.
%
%   Output:
%       height_mm       : Required vessel height (length) [mm] at each time
%                         point so that the cylindrical volume equals the
%                         provided contrast volume.
%
%   Notes:
%       • All inputs must be vectors of identical length.
%       • A strictly positive radius is required to avoid division by zero.
%
%   Example:
%       t_s = linspace(0, 5, 6);
%       radius_mm = 2.5 * ones(size(t_s));
%       V_contrast_mm3 = [0 10 25 50 70 90];
%       height_mm = computeContrastWashIn(t_s, radius_mm, V_contrast_mm3);
%       % height_mm gives the necessary enhanced length at each time point.

arguments
    t_s (:) double {mustBeFinite}
    radius_mm (:) double {mustBePositive, mustBeFinite}
    V_contrast_mm3 (:) double {mustBeFinite, mustBeNonnegative}
end

numSamples = numel(t_s);
if numel(radius_mm) ~= numSamples || numel(V_contrast_mm3) ~= numSamples
    error('computeContrastWashIn:SizeMismatch', ...
        't_s, radius_mm, and V_contrast_mm3 must have identical lengths.');
end

% Cylinder volume: V = pi * r^2 * h  =>  h = V / (pi * r^2)
area_mm2 = pi .* (radius_mm .^ 2);
height_mm = V_contrast_mm3 ./ area_mm2;
end
