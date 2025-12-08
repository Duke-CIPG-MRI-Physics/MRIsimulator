classdef AnalyticalEllipticalCylinder3D < AnalyticalShape3D
    % AnalyticalEllipticalCylinder3D
    %   Finite cylinder with elliptical cross-section and analytic 3D FT.
    %
    %   Body-frame definition:
    %     - Elliptical cross-section: x^2/a^2 + y^2/b^2 <= 1
    %     - Axis along z
    %     - Length L_mm (z ∈ [-L_mm/2, +L_mm/2])
    %
    %   World-frame pose handled by AnalyticalShape3D:
    %     - center: [x0,y0,z0] in mm
    %     - rollPitchYaw: [roll pitch yaw] in degrees, body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)

    %% Private geometry parameters
    properties (Access = private)
        a_mm = 5;    % semi-axis in x (mm), numeric or time-varying spec
        b_mm = 3;    % semi-axis in y (mm), numeric or time-varying spec
        L_mm = 20;   % length in z (mm), numeric or time-varying spec
    end

    %% Constructor
    methods
        function obj = AnalyticalEllipticalCylinder3D(a_mm, b_mm, L_mm, intensity, center, rollPitchYaw)
            % Constructor
            %
            %   obj = AnalyticalEllipticalCylinder3D()
            %   obj = AnalyticalEllipticalCylinder3D(a_mm, b_mm, L_mm)
            %   obj = AnalyticalEllipticalCylinder3D(a_mm, b_mm, L_mm, intensity)
            %   obj = AnalyticalEllipticalCylinder3D(a_mm, b_mm, L_mm, intensity, center)
            %   obj = AnalyticalEllipticalCylinder3D(a_mm, b_mm, L_mm, intensity, center, rollPitchYaw)

            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);

            if nargin >= 1 && ~isempty(a_mm)
                obj.setA(a_mm);
            end
            if nargin >= 2 && ~isempty(b_mm)
                obj.setB(b_mm);
            end
            if nargin >= 3 && ~isempty(L_mm)
                obj.setLength(L_mm);
            end
        end
    end

    %% Public geometry getters/setters
    methods
        function a = getA(obj)
            a = obj.evaluateParameter(obj.a_mm, 'a');
        end

        function setA(obj, newA)
            % setA  Accepts numeric values or @(t)->a_mm waveforms.
            arguments
                obj
                newA
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            aSpec = obj.normalizeGeometryInput(newA, validator, 'a');

            if ~isequal(obj.a_mm, aSpec)
                obj.a_mm = aSpec;
                obj.markShapeChanged();
            end
        end

        function b = getB(obj)
            b = obj.evaluateParameter(obj.b_mm, 'b');
        end

        function setB(obj, newB)
            % setB  Accepts numeric values or @(t)->b_mm waveforms.
            arguments
                obj
                newB
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            bSpec = obj.normalizeGeometryInput(newB, validator, 'b');

            if ~isequal(obj.b_mm, bSpec)
                obj.b_mm = bSpec;
                obj.markShapeChanged();
            end
        end

        function L = getLength(obj)
            L = obj.evaluateParameter(obj.L_mm, 'L');
        end

        function setLength(obj, newL)
            % setLength  Accepts numeric values or @(t)->L_mm waveforms.
            arguments
                obj
                newL
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            lSpec = obj.normalizeGeometryInput(newL, validator, 'L');

            if ~isequal(obj.L_mm, lSpec)
                obj.L_mm = lSpec;
                obj.markShapeChanged();
            end
        end

    end

    %% Analytic BODY-frame FT and inside-test
    methods (Access = protected)
        function S = kspaceBodyGeometry(obj, kx_body, ky_body, kz_body)
            % kspaceBodyGeometry
            %   Analytic FT of elliptical cylinder aligned with body z-axis.
            %
            %   Using affine-scaling relationship from circular cylinder:
            %     - A = diag(a, b, 1)
            %     - det(A) = a*b
            %     - k' = A^T k => k_perp' = sqrt((a kx)^2 + (b ky)^2)
            %
            %   Circular cylinder FT in terms of k_perp':
            %     radial(k_perp') = R * J1(2*pi*R*k_perp') / k_perp'   (R = 1 here)
            %                      -> pi * R^2 at k_perp' = 0
            %
            %   Here we set base radius R0 = 1 and scale volume by a*b
            %   to get the elliptical cross-section area = pi*a*b.

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('AnalyticalEllipticalCylinder3D:kspaceBodyGeometry:SizeMismatch', ...
                    'kx_body, ky_body, kz_body must have identical sizes.');
            end

            a = obj.requireScalarOrSize(obj.getA(), kx_body, 'a');
            b = obj.requireScalarOrSize(obj.getB(), kx_body, 'b');
            L = obj.requireScalarOrSize(obj.getLength(), kx_body, 'L');

            % Scaled radial frequency
            k_perp_prime = sqrt( (a .* kx_body).^2 + (b .* ky_body).^2 );

            % Radial part (base cylinder radius R0 = 1, then scale by a*b)
            radial = zeros(size(k_perp_prime));
            idx_r = (k_perp_prime ~= 0);
            if any(idx_r(:))
                x_r = 2*pi .* k_perp_prime(idx_r);    % R0=1
                radial(idx_r) = (a(idx_r) .* b(idx_r)) .* besselj(1, x_r) ./ k_perp_prime(idx_r);
            end
            % k_perp' -> 0 limit: area of ellipse = pi*a*b
            radial(~idx_r) = pi .* a(~idx_r) .* b(~idx_r);

            % Axial part (same as finite cylinder of length L)
            kz = kz_body;
            zfac = zeros(size(kz));
            idx_z = (kz ~= 0);
            if any(idx_z(:))
                x_z = pi .* L(idx_z) .* kz(idx_z);
                zfac(idx_z) = L(idx_z) .* (sin(x_z) ./ x_z);   % L * sinc(pi*L*kz)
            end
            % kz -> 0 limit: length
            zfac(~idx_z) = L(~idx_z);

            S = radial .* zfac;
        end

        function frac = percentInsideShape(obj, xb, yb, zb)
            % percentInsideShape
            %   BODY-frame inside test for elliptical cylinder.

            a = obj.requireScalarOrSize(obj.getA(), xb, 'a');
            b = obj.requireScalarOrSize(obj.getB(), xb, 'b');
            L = obj.requireScalarOrSize(obj.getLength(), xb, 'L');

            inside = (xb./a).^2 + (yb./b).^2 <= 1 & (abs(zb) <= L./2);
            frac = double(inside);
        end
    end
end
