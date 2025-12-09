classdef AnalyticalEllipticalCylinder3D < AnalyticalShape3D
    % AnalyticalEllipticalCylinder3D
    %   Finite cylinder with an elliptical cross section aligned to BODY
    %   axes.

    properties (Access = protected)
        a_mm double {mustBeNonnegative} = 1; % semi-axis along x
        b_mm double {mustBeNonnegative} = 1; % semi-axis along y
        length_mm double {mustBeNonnegative} = 1;
    end

    methods
        function obj = AnalyticalEllipticalCylinder3D(a_mm, b_mm, length_mm, intensity, center, rollPitchYaw)
            if nargin < 6
                rollPitchYaw = [0 0 0];
            end
            if nargin < 5
                center = [0 0 0];
            end
            if nargin < 4
                intensity = 1;
            end

            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);
            if nargin >= 1 && ~isempty(a_mm)
                obj.a_mm = a_mm;
            end
            if nargin >= 2 && ~isempty(b_mm)
                obj.b_mm = b_mm;
            end
            if nargin >= 3 && ~isempty(length_mm)
                obj.length_mm = length_mm;
            end
        end

        function setAxes(obj, a_mm, b_mm)
            obj.a_mm = a_mm;
            obj.b_mm = b_mm;
        end

        function setLength(obj, length_mm)
            obj.length_mm = length_mm;
        end

        function vol = calculateVolume(obj)
            vol = pi .* obj.a_mm .* obj.b_mm .* obj.length_mm;
        end
    end

    methods (Access = protected)
        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            kr = sqrt((obj.a_mm .* kx).^2 + (obj.b_mm .* ky).^2);
            arg = 2 * pi .* kr;

            radial = (obj.a_mm .* obj.b_mm .* besselj(1, arg)) ./ kr;
            radial(kr == 0) = pi .* obj.a_mm .* obj.b_mm;

            axial = obj.length_mm .* sinc(obj.length_mm .* kz);

            S_body = radial .* axial;
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            inEllipse = (xb ./ obj.a_mm).^2 + (yb ./ obj.b_mm).^2 <= 1;
            inHeight = abs(zb) <= obj.length_mm ./ 2;
            percent = double(inEllipse & inHeight);
        end
    end
end
