classdef AnalyticalBox3D < AnalyticalShape3D
    % AnalyticalBox3D
    %   Solid axis-aligned rectangular box with analytic 3D FT.
    %
    %   BODY frame definition:
    %     - Centered at origin
    %     - x ∈ [-Lx/2, +Lx/2]
    %     - y ∈ [-Ly/2, +Ly/2]
    %     - z ∈ [-Lz/2, +Lz/2]
    %
    %   WORLD pose handled by AnalyticalShape3D:
    %     - center: [x0,y0,z0] in mm
    %     - rollPitchYaw: [roll pitch yaw] in degrees, body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)

    %% Private geometry parameters
    properties (Access = private)
        Lx_mm = 10;   % box length in x (mm)
        Ly_mm = 10;   % box length in y (mm)
        Lz_mm = 10;   % box length in z (mm)
    end

    %% Constructor
    methods
        function obj = AnalyticalBox3D(Lx_mm, Ly_mm, Lz_mm, intensity, center, rollPitchYaw)
            % Constructor
            %
            %   obj = AnalyticalBox3D()
            %   obj = AnalyticalBox3D(Lx, Ly, Lz)
            %   obj = AnalyticalBox3D(Lx, Ly, Lz, intensity)
            %   obj = AnalyticalBox3D(Lx, Ly, Lz, intensity, center)
            %   obj = AnalyticalBox3D(Lx, Ly, Lz, intensity, center, rollPitchYaw)

            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);

            if nargin >= 1 && ~isempty(Lx_mm)
                obj.setLengthX(Lx_mm);
            end
            if nargin >= 2 && ~isempty(Ly_mm)
                obj.setLengthY(Ly_mm);
            end
            if nargin >= 3 && ~isempty(Lz_mm)
                obj.setLengthZ(Lz_mm);
            end
        end
    end

    %% Public geometry getters/setters
    methods
        function Lx = getLengthX(obj)
            Lx = obj.evaluateParameter(obj.Lx_mm, 'Lx');
        end

        function setLengthX(obj, newLx, opts)
            arguments
                obj
                newLx
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            LxSpec = obj.normalizeGeometryInput(newLx, validator, 'Lx');

            if ~isequal(obj.Lx_mm, LxSpec)
                obj.Lx_mm = LxSpec;
                obj.markShapeChanged();
            end
        end

        function Ly = getLengthY(obj)
            Ly = obj.evaluateParameter(obj.Ly_mm, 'Ly');
        end

        function setLengthY(obj, newLy, opts)
            arguments
                obj
                newLy
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            LySpec = obj.normalizeGeometryInput(newLy, validator, 'Ly');

            if ~isequal(obj.Ly_mm, LySpec)
                obj.Ly_mm = LySpec;
                obj.markShapeChanged();
            end
        end

        function Lz = getLengthZ(obj)
            Lz = obj.evaluateParameter(obj.Lz_mm, 'Lz');
        end

        function setLengthZ(obj, newLz, opts)
            arguments
                obj
                newLz
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            LzSpec = obj.normalizeGeometryInput(newLz, validator, 'Lz');

            if ~isequal(obj.Lz_mm, LzSpec)
                obj.Lz_mm = LzSpec;
                obj.markShapeChanged();
            end
        end

    end

    %% Analytic BODY-frame FT and inside-test
    methods (Access = protected)
        function S = kspaceBodyGeometry(obj, kx_body, ky_body, kz_body)
            % kspaceBodyGeometry
            %   Analytic FT of a solid box centered at origin in body frame.
            %
            %   Uses MATLAB's sinc(u) = sin(pi*u)/(pi*u).
            %   For 1D box of length L:
            %       F(k) = L * sinc(k*L)
            %
            %   So in 3D:
            %       F(kx,ky,kz) = Lx*Ly*Lz *
            %                     sinc(kx*Lx) * sinc(ky*Ly) * sinc(kz*Lz)

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('AnalyticalBox3D:kspaceBodyGeometry:SizeMismatch', ...
                    'kx_body, ky_body, kz_body must have identical sizes.');
            end

            Lx = obj.requireScalarOrSize(obj.getLengthX(), kx_body, 'Lx');
            Ly = obj.requireScalarOrSize(obj.getLengthY(), kx_body, 'Ly');
            Lz = obj.requireScalarOrSize(obj.getLengthZ(), kx_body, 'Lz');

            Sx = Lx .* sinc(kx_body .* Lx);
            Sy = Ly .* sinc(ky_body .* Ly);
            Sz = Lz .* sinc(kz_body .* Lz);

            S = Sx .* Sy .* Sz;
        end

        function frac = percentInsideShape(obj, xb, yb, zb)
            % percentInsideShape
            %   BODY-frame inside test for rectangular box.

            Lx = obj.requireScalarOrSize(obj.getLengthX(), xb, 'Lx');
            Ly = obj.requireScalarOrSize(obj.getLengthY(), xb, 'Ly');
            Lz = obj.requireScalarOrSize(obj.getLengthZ(), xb, 'Lz');

            inside = (abs(xb) <= Lx./2) & (abs(yb) <= Ly./2) & (abs(zb) <= Lz./2);
            frac = double(inside);
        end
    end
end
