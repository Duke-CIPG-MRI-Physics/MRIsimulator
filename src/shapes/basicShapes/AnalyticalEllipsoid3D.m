classdef AnalyticalEllipsoid3D < AnalyticalShape3D
    % AnalyticalEllipsoid3D
    %   Solid ellipsoid aligned with BODY axes.

    properties (Access = protected)
        a_mm double {mustBeNonnegative} = 1;
        b_mm double {mustBeNonnegative} = 1;
        c_mm double {mustBeNonnegative} = 1;
    end

    methods
        function obj = AnalyticalEllipsoid3D(a_mm, b_mm, c_mm, intensity, center, rollPitchYaw)
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
            if nargin >= 3 && ~isempty(c_mm)
                obj.c_mm = c_mm;
            end
        end

        function setAxes(obj, a_mm, b_mm, c_mm)
            obj.a_mm = a_mm;
            obj.b_mm = b_mm;
            obj.c_mm = c_mm;
        end

        function vol = calculateVolume(obj)
            vol = (4/3) * pi .* obj.a_mm .* obj.b_mm .* obj.c_mm;
        end
    end

    methods (Access = protected)
        function S_body = bodyKspace(obj, kx, ky, kz)
            kScaled = sqrt((obj.a_mm .* kx).^2 + (obj.b_mm .* ky).^2 + (obj.c_mm .* kz).^2);
            arg = 2 * pi .* kScaled;

            % Spherical Bessel j1(x) = (sin(x) - x cos(x)) / x^2
            numerator = sin(arg) - arg .* cos(arg);
            denom = (arg).^3;

            S_body = obj.calculateVolume() .* 3 .* numerator ./ denom;
            S_body(kScaled == 0) = obj.calculateVolume();
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            scaled = (xb ./ obj.a_mm).^2 + (yb ./ obj.b_mm).^2 + (zb ./ obj.c_mm).^2;
            percent = double(scaled <= 1);
        end
    end
end
