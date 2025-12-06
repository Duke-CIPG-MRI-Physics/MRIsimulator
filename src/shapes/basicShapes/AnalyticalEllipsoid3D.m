classdef AnalyticalEllipsoid3D < AnalyticalShape3D
    % AnalyticalEllipsoid3D
    %   Solid ellipsoid with analytic 3D FT in k-space.
    %
    %   BODY frame definition:
    %     (x/a_mm)^2 + (y/b_mm)^2 + (z/c_mm)^2 <= 1
    %
    %   WORLD pose handled by AnalyticalShape3D:
    %     - center: [x0,y0,z0] in mm
    %     - rollPitchYaw: [roll pitch yaw] in degrees, body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)

    %% Private geometry parameters
    properties (Access = private)
        a_mm = 5;   % semi-axis in x (mm), numeric or time-varying spec
        b_mm = 5;   % semi-axis in y (mm), numeric or time-varying spec
        c_mm = 5;   % semi-axis in z (mm), numeric or time-varying spec
    end

    %% Constructor
    methods
        function obj = AnalyticalEllipsoid3D(a_mm, b_mm, c_mm, intensity, center, rollPitchYaw)
            % Constructor
            %
            %   obj = AnalyticalEllipsoid3D()
            %   obj = AnalyticalEllipsoid3D(a_mm, b_mm, c_mm)
            %   obj = AnalyticalEllipsoid3D(a_mm, b_mm, c_mm, intensity)
            %   obj = AnalyticalEllipsoid3D(a_mm, b_mm, c_mm, intensity, center)
            %   obj = AnalyticalEllipsoid3D(a_mm, b_mm, c_mm, intensity, center, rollPitchYaw)

            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);

            if nargin >= 1 && ~isempty(a_mm)
                obj.setA(a_mm);
            end
            if nargin >= 2 && ~isempty(b_mm)
                obj.setB(b_mm);
            end
            if nargin >= 3 && ~isempty(c_mm)
                obj.setC(c_mm);
            end
        end
    end

    %% Public geometry getters/setters
    methods
        function a = getA(obj)
            a = obj.a_mm;
        end

        function setA(obj, newA, opts)
            % setA  Accepts numeric values or @(t)->a_mm waveforms.
            arguments
                obj
                newA
                opts.Cache logical = true
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            aSpec = obj.normalizeGeometryInput(newA, validator, opts.Cache, 'a');

            if ~isequal(obj.a_mm, aSpec)
                obj.a_mm = aSpec;
                obj.markShapeChanged();
            end
        end

        function b = getB(obj)
            b = obj.b_mm;
        end

        function setB(obj, newB, opts)
            % setB  Accepts numeric values or @(t)->b_mm waveforms.
            arguments
                obj
                newB
                opts.Cache logical = true
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            bSpec = obj.normalizeGeometryInput(newB, validator, opts.Cache, 'b');

            if ~isequal(obj.b_mm, bSpec)
                obj.b_mm = bSpec;
                obj.markShapeChanged();
            end
        end

        function c = getC(obj)
            c = obj.c_mm;
        end

        function setC(obj, newC, opts)
            % setC  Accepts numeric values or @(t)->c_mm waveforms.
            arguments
                obj
                newC
                opts.Cache logical = true
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            cSpec = obj.normalizeGeometryInput(newC, validator, opts.Cache, 'c');

            if ~isequal(obj.c_mm, cSpec)
                obj.c_mm = cSpec;
                obj.markShapeChanged();
            end
        end

    end

    %% Analytic BODY-frame FT and inside-test
    methods (Access = protected)
        function S = kspaceBodyGeometry(obj, kx_body, ky_body, kz_body)
            % kspaceBodyGeometry
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
                error('AnalyticalEllipsoid3D:kspaceBodyGeometry:SizeMismatch', ...
                    'kx_body, ky_body, kz_body must have identical sizes.');
            end

            a = obj.requireScalarOrSize(obj.a_mm, kx_body, 'a');
            b = obj.requireScalarOrSize(obj.b_mm, kx_body, 'b');
            c = obj.requireScalarOrSize(obj.c_mm, kx_body, 'c');

            % "Effective" radial frequency in scaled coordinates
            kprime = sqrt((a .* kx_body).^2 + ...
                          (b .* ky_body).^2 + ...
                          (c .* kz_body).^2);

            S = zeros(size(kprime));

            idx = (kprime ~= 0);
            volume = a .* b .* c;

            if any(idx(:))
                x = 2*pi .* kprime(idx);   % dimensionless
                % Base sphere (R=1) FT
                Fs = 4*pi .* (sin(x) - x .* cos(x)) ./ (x.^3);
                % Scale by volume factor a*b*c
                S(idx) = volume(idx) .* Fs;
            end

            % k' -> 0 limit: ellipsoid volume
            S(~idx) = (4/3) * pi .* volume(~idx);
        end

        function frac = percentInsideShape(obj, xb, yb, zb)
            % percentInsideShape
            %   BODY-frame inside test for ellipsoid.

            a = obj.requireScalarOrSize(obj.a_mm, xb, 'a');
            b = obj.requireScalarOrSize(obj.b_mm, xb, 'b');
            c = obj.requireScalarOrSize(obj.c_mm, xb, 'c');

            inside = (xb./a).^2 + (yb./b).^2 + (zb./c).^2 <= 1;
            frac = double(inside);
        end
    end
end
