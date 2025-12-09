classdef AnalyticalSphere3D < AnalyticalEllipsoid3D
    % AnalyticalSphere3D
    %   Convenience wrapper around AnalyticalEllipsoid3D with equal axes.

    methods
        function obj = AnalyticalSphere3D(radius_mm, intensity, center, rollPitchYaw)
            if nargin < 4
                rollPitchYaw = [0 0 0];
            end
            if nargin < 3
                center = [0 0 0];
            end
            if nargin < 2
                intensity = 1;
            end

            obj@AnalyticalEllipsoid3D(radius_mm, radius_mm, radius_mm, ...
                intensity, center, rollPitchYaw);
        end
    end
end
