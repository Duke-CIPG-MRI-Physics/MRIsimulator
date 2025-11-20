classdef AnalyticalEllipsoid3D < AnalyticalShape3D
    % AnalyticalEllipsoid3D
    %   Solid ellipsoid with analytic 3D FT in k-space.
    %
    %   BODY frame definition:
    %     (x/a_mm)^2 + (y/b_mm)^2 + (z/c_mm)^2 <= 1
    %
    %   WORLD pose handled by AnalyticalShape3D:
    %     - center: [x0,y0,z0] in mm
    %     - roll_deg, pitch_deg, yaw_deg: orientation (deg), body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)

    %% Private geometry parameters
    properties (Access = private)
        a_mm (1,1) double {mustBePositive} = 5;   % semi-axis in x (mm)
        b_mm (1,1) double {mustBePositive} = 5;   % semi-axis in y (mm)
        c_mm (1,1) double {mustBePositive} = 5;   % semi-axis in z (mm)
    end

    %% Constructor
    methods
        function obj = AnalyticalEllipsoid3D(center, roll_deg, pitch_deg, yaw_deg, a_mm, b_mm, c_mm)
            % Constructor
            %
            %   obj = AnalyticalEllipsoid3D()
            %   obj = AnalyticalEllipsoid3D(center)
            %   obj = AnalyticalEllipsoid3D(center, roll, pitch, yaw)
            %   obj = AnalyticalEllipsoid3D(center, roll, pitch, yaw, a_mm, b_mm, c_mm)

            obj@AnalyticalShape3D();

            if nargin >= 1 && ~isempty(center)
                obj.setCenter(center);
            end

            if nargin >= 4 && ~isempty(roll_deg)
                obj.setOrientation(roll_deg, pitch_deg, yaw_deg);
            end

            if nargin >= 5 && ~isempty(a_mm)
                obj.setA(a_mm);
            end
            if nargin >= 6 && ~isempty(b_mm)
                obj.setB(b_mm);
            end
            if nargin >= 7 && ~isempty(c_mm)
                obj.setC(c_mm);
            end
        end
    end

    %% Public geometry getters/setters
    methods
        function a = getA(obj)
            a = obj.a_mm;
        end

        function setA(obj, newA)
            arguments
                obj
                newA (1,1) double {mustBePositive}
            end

            if obj.a_mm ~= newA
                obj.a_mm = newA;
                obj.markShapeChanged();
            end
        end

        function b = getB(obj)
            b = obj.b_mm;
        end

        function setB(obj, newB)
            arguments
                obj
                newB (1,1) double {mustBePositive}
            end

            if obj.b_mm ~= newB
                obj.b_mm = newB;
                obj.markShapeChanged();
            end
        end

        function c = getC(obj)
            c = obj.c_mm;
        end

        function setC(obj, newC)
            arguments
                obj
                newC (1,1) double {mustBePositive}
            end

            if obj.c_mm ~= newC
                obj.c_mm = newC;
                obj.markShapeChanged();
            end
        end

        function volume_mm3 = calculateVolume(obj)
            a = obj.a_mm;
            b = obj.b_mm;
            c = obj.c_mm;
            volume_mm3 = (4/3) * pi * a * b * c;
        end
    end

    %% Analytic BODY-frame FT and inside-test
    methods (Access = protected)
        function S = bodyKspace(obj, kx_body, ky_body, kz_body)
            % bodyKspace
            %   Analytic FT of a solid ellipsoid centered at origin in body frame.
            %
            %   Uses affine scaling from a unit sphere:
            %     k' = sqrt((a*kx)^2 + (b*ky)^2 + (c*kz)^2)
            %     F_ellip(k) = a*b*c * F_sphere_R1(k')
            %
            %   where
            %     F_sphere_R1(k') = 4*pi * [sin(2*pi*k') - 2*pi*k' * cos(2*pi*k')]
            %                        / (2*pi*k')^3
            %
            %   At k' = 0, F_ellip = 4/3*pi*a*b*c (ellipsoid volume).

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('AnalyticalEllipsoid3D:bodyKspace:SizeMismatch', ...
                    'kx_body, ky_body, kz_body must have identical sizes.');
            end

            a = obj.a_mm;
            b = obj.b_mm;
            c = obj.c_mm;

            % "Effective" radial frequency in scaled coordinates
            kprime = sqrt((a .* kx_body).^2 + ...
                          (b .* ky_body).^2 + ...
                          (c .* kz_body).^2);

            S = zeros(size(kprime));

            idx = (kprime ~= 0);
            if any(idx(:))
                x = 2*pi .* kprime(idx);   % dimensionless
                % Base sphere (R=1) FT
                Fs = 4*pi .* (sin(x) - x .* cos(x)) ./ (x.^3);
                % Scale by volume factor a*b*c
                S(idx) = (a * b * c) .* Fs;
            end

            % k' -> 0 limit: ellipsoid volume
            S(~idx) = 4/3 * pi * a * b * c;
        end

        function frac = percentInsideShape(obj, xb, yb, zb)
            % percentInsideShape
            %   BODY-frame inside test for ellipsoid.

            a = obj.a_mm;
            b = obj.b_mm;
            c = obj.c_mm;

            inside = (xb./a).^2 + (yb./b).^2 + (zb./c).^2 <= 1;
            frac = double(inside);
        end
    end
end
