classdef AnalyticalCylinder3D < AnalyticalShape3D
    % AnalyticalCylinder3D
    %   Solid finite cylinder with analytic 3D Fourier transform (BODY frame).
    %
    %   BODY frame definition:
    %       - Cylinder axis aligned with +z
    %       - Radius = R_mm        (mm)
    %       - Length = L_mm        (mm), z ∈ [-L_mm/2, +L_mm/2]
    %
    %   WORLD pose:
    %       - center via translation
    %       - orientation applied via rotation matrix from AnalyticalShape3D
    %
    %   BODY-frame FT:
    %       k_perp = sqrt(kx^2 + ky^2)
    %
    %       radial(k_perp) =
    %         R * J1(2πR k_perp) / k_perp        (k_perp ≠ 0)
    %         π R^2                              (k_perp = 0)
    %
    %       zfac(kz) =
    %         L * sin(πL kz) / (πL kz)           (kz ≠ 0)
    %         L                                  (kz = 0)
    %
    %       S(k) = radial(k_perp) .* zfac(kz)
    %
    %   Image-domain helpers:
    %       estimateImageShape(xMesh,yMesh,zMesh)

    %% Private geometry parameters
    properties (Access = private)
        R_mm double {mustBePositive} = 5;   % radius [mm]
        L_mm double {mustBeNonnegative} = 20;  % length [mm]
    end

    %% Constructor
    methods
        function obj = AnalyticalCylinder3D(R_mm, L_mm, intensity, center, rollPitchYaw)
            % Constructor follows the same pattern as AnalyticalSphere3D
            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);  % must be first, unconditional

            % Optional geometry
            if nargin >= 1 && ~isempty(R_mm)
                obj.setRadius(R_mm);
            end
            if nargin >= 2 && ~isempty(L_mm)
                obj.setLength(L_mm);
            end
        end
    end

    %% Public geometry getters/setters (fires shapeChanged)
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

        function L = getLength(obj)
            L = obj.L_mm;
        end

        function setLength(obj, newLength)
            arguments
                obj
                newLength double {mustBeNonnegative}
            end
            if ~isequal(obj.L_mm, newLength)
                obj.L_mm = newLength;
                obj.markShapeChanged();
            end
        end

    end

    %% BODY-frame analytic FT and inside-test
    methods (Access = protected)

        % ----------- Analytic Fourier transform (BODY frame) -------------
        function S = kspaceBodyGeometry(obj, kx_body, ky_body, kz_body)
            arguments
                obj
                kx_body double
                ky_body double
                kz_body double
            end

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('AnalyticalCylinder3D:SizeMismatch');
            end

            R = obj.requireScalarOrSize(obj.R_mm, kx_body, 'R');
            L = obj.requireScalarOrSize(obj.L_mm, kx_body, 'L');

            % Degenerate cylinder (zero radius or length) has zero volume → zero-valued FT
            if all(R(:) == 0 | L(:) == 0)
                S = zeros(size(kx_body));
                return;
            end

            % Radial term
            k_perp = sqrt(kx_body.^2 + ky_body.^2);
            radial = zeros(size(k_perp));

            nonZero_idx = (k_perp ~= 0);
            if any(nonZero_idx(:))
                x_r = 2*pi .* R(nonZero_idx) .* k_perp(nonZero_idx);
                radial(nonZero_idx) = R(nonZero_idx) .* besselj(1, x_r) ./ k_perp(nonZero_idx);
            end
            radial(~nonZero_idx) = pi .* (R(~nonZero_idx).^2);   % limit as k_perp → 0

            % Axial term
            zfac = zeros(size(kz_body));
            nonZero_idx = (kz_body ~= 0) & (L ~= 0);
            if any(nonZero_idx(:))
                x_z = pi .* L(nonZero_idx) .* kz_body(nonZero_idx);
                zfac(nonZero_idx) = L(nonZero_idx) .* (sin(x_z) ./ x_z);
            end
            zfac(~nonZero_idx) = L(~nonZero_idx);

            S = radial .* zfac;
        end

        % ----------------- Inside test (BODY frame) ----------------------
        function frac = percentInsideShape(obj, xb, yb, zb)
            % 0/1 mask (upgradeable later to partial volume)
            R = obj.requireScalarOrSize(obj.R_mm, xb, 'R');
            L = obj.requireScalarOrSize(obj.L_mm, xb, 'L');

            inside = (xb.^2 + yb.^2 <= R.^2) & (abs(zb) <= L./2);
            frac = double(inside);
        end
    end
end
