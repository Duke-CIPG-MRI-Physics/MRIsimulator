classdef AnalyticalSphere3D < AnalyticalShape3D
    % AnalyticalSphere3D
    %   Solid sphere with analytic 3D Fourier transform in the BODY frame.
    %
    %   BODY frame:
    %       - Sphere centered at (0,0,0)
    %       - Radius = R_mm   (mm)
    %
    %   WORLD pose:
    %       - center is applied (translation)
    %       - orientation is IGNORED (sphere is isotropic)
    %
    %   FT:
    %       k = sqrt(kx^2 + ky^2 + kz^2), x = 2π k R
    %       S(k) = 4π R^3 [ sin(x) - x cos(x) ] / x^3
    %       S(0) = 4/3 π R^3
    %
    %   Image-domain helper:
    %       estimateImageShape(xMesh,yMesh,zMesh)   → fraction map

    %% Private geometry parameter
    properties (Access = private)
        R_mm double {mustBePositive} = 5;   % sphere radius [mm]
    end

    %% Constructor
    methods
        function obj = AnalyticalSphere3D(R_mm, intensity, center, rollPitchYaw)
            % Constructor
            %
            %   obj = AnalyticalSphere3D()
            %   obj = AnalyticalSphere3D(R_mm)
            %   obj = AnalyticalSphere3D(R_mm, intensity)
            %   obj = AnalyticalSphere3D(R_mm, intensity, center)
            %   obj = AnalyticalSphere3D(R_mm, intensity, center, rollPitchYaw)

            % Base construction (pose + intensity)
            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);

            % Radius (geometry matters)
            if nargin >= 1 && ~isempty(R_mm)
                obj.setRadius(R_mm);
            end
        end
    end

    %% Override orientation: ignore any requested rotation
    methods
        function setOrientation(obj, ~, ~, ~)
            % For a sphere, orientation does not change the shape or its FT.
            obj.setRollPitchYaw([0 0 0]);
            % No markShapeChanged(): orientation irrelevant for sphere.
        end
    end

    %% Public geometry API
    methods
        function R = getRadius(obj)
            R = obj.R_mm;
        end

        function setRadius(obj, newRadius)
            arguments
                obj
                newRadius double {mustBePositive}
            end

            if ~isequal(obj.R_mm, newRadius)
                obj.R_mm = newRadius;
                obj.markShapeChanged();
            end
        end

    end

    %% Analytic BODY-frame FT
    methods (Access = protected)
        function S = kspaceBodyGeometry(obj, kx_body, ky_body, kz_body)
            arguments
                obj
                kx_body double
                ky_body double
                kz_body double
            end

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('AnalyticalSphere3D:SizeMismatch', ...
                    'kx_body, ky_body, kz_body must be same size.');
            end

            k = sqrt(kx_body.^2 + ky_body.^2 + kz_body.^2);
            R = obj.requireScalarOrSize(obj.R_mm, kx_body, 'R');
            x = 2*pi .* R .* k;

            S = zeros(size(k));

            idx = (k ~= 0);
            if any(idx(:))
                x_nz = x(idx);
                R_nz = R(idx);
                S(idx) = 4*pi .* (R_nz.^3) .* (sin(x_nz) - x_nz .* cos(x_nz)) ./ (x_nz.^3);
            end

            % k = 0 → sphere volume
            volume = (4/3) * pi .* (R.^3);
            S(~idx) = volume(~idx);
        end

        function imageShape = percentInsideShape(obj, xb, yb, zb)
            % percentInsideShape (sphere)
            %   BODY-frame fraction of each pixel inside the sphere.
            %
            %   Current implementation:
            %       - 0/1 mask: 1 if pixel center inside sphere, 0 otherwise.
            %       - This is a placeholder; you can upgrade to partial volume
            %         by supersampling in BODY coordinates.
            %
            %   Inputs:
            %       xb,yb,zb : BODY coords of pixel centers [mm]
            %       varargin : reserved for options (e.g., nSubSamples)
            %
            %   Output:
            %       imageShape : double in [0,1], same size as xb.

            %#ok<*INUSD>  % suppress unused varargin warning for now

            R = obj.requireScalarOrSize(obj.R_mm, xb, 'R');
            inside = (xb.^2 + yb.^2 + zb.^2) <= R.^2;

            imageShape = double(inside);   % 0 or 1 (no partial volume yet)
        end
    end
end
