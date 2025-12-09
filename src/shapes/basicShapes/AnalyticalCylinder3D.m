classdef AnalyticalCylinder3D < AnalyticalShape3D
    % AnalyticalCylinder3D
    %   Finite cylinder aligned with the BODY +z axis.

    properties (Access = protected)
        radius_mm double {mustBeNonnegative} = 1;
        length_mm double {mustBeNonnegative} = 1;
    end

    methods
        function obj = AnalyticalCylinder3D(radius_mm, length_mm, intensity, center, rollPitchYaw)
            if nargin < 5
                rollPitchYaw = [0 0 0];
            end
            if nargin < 4
                center = [0 0 0];
            end
            if nargin < 3
                intensity = 1;
            end

            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);
            if nargin >= 1 && ~isempty(radius_mm)
                obj.radius_mm = radius_mm;
            end
            if nargin >= 2 && ~isempty(length_mm)
                obj.length_mm = length_mm;
            end
        end

        function setRadius(obj, radius_mm)
            obj.radius_mm = radius_mm;
            obj.markShapeChanged();
        end

        function setLength(obj, length_mm)
            obj.length_mm = length_mm;
            obj.markShapeChanged();
        end

        function r = getRadius(obj)
            r = obj.radius_mm;
        end

        function L = getLength(obj)
            L = obj.length_mm;
        end

        function vol = calculateVolume(obj)
            vol = pi .* (obj.radius_mm.^2) .* obj.length_mm;
        end
    end

    methods (Access = protected)
        function S_body = bodyKspace(obj, kx, ky, kz)
            kr = sqrt(kx.^2 + ky.^2);
            arg = 2*pi*obj.radius_mm .* kr;

            radial = (obj.radius_mm .* besselj(1, arg)) ./ kr;
            radial(kr == 0) = pi .* obj.radius_mm.^2;

            axial = obj.length_mm .* sinc(obj.length_mm .* kz);

            S_body = radial .* axial;
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            inRadius = (xb.^2 + yb.^2) <= (obj.radius_mm.^2);
            inHeight = abs(zb) <= (obj.length_mm ./ 2);
            percent = double(inRadius & inHeight);
        end
    end
end
