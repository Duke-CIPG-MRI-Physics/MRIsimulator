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
    %       calculateImageOutline(xMesh,yMesh,zMesh)

    %% Private geometry parameters
    properties (Access = private)
        R_mm (1,1) double {mustBePositive} = 5;   % radius [mm]
        L_mm (1,1) double {mustBePositive} = 20;  % length [mm]
    end

    %% Constructor
    methods
        function obj = AnalyticalCylinder3D(center, roll_deg, pitch_deg, yaw_deg, R_mm, L_mm)
            % Constructor follows the same pattern as AnalyticalSphere3D
            obj@AnalyticalShape3D();  % must be first, unconditional

            % Optional center
            if nargin >= 1 && ~isempty(center)
                obj.setCenter(center);
            end

            % Optional orientation
            if nargin >= 4 && ~isempty(roll_deg)
                obj.setOrientation(roll_deg, pitch_deg, yaw_deg);
            end

            % Optional geometry
            if nargin >= 5 && ~isempty(R_mm)
                obj.setRadius(R_mm);
            end
            if nargin >= 6 && ~isempty(L_mm)
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
                newRadius (1,1) double {mustBePositive}
            end
            if obj.R_mm ~= newRadius
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
                newLength (1,1) double {mustBePositive}
            end
            if obj.L_mm ~= newLength
                obj.L_mm = newLength;
                obj.markShapeChanged();
            end
        end
    end

    %% Outline extraction (same design as sphere)
    methods
        function [xOut, yOut, zOut, insideMask] = calculateImageOutline(obj, xMesh, yMesh, zMesh)
            arguments
                obj
                xMesh double
                yMesh double
                zMesh double
            end

            if ~isequal(size(xMesh), size(yMesh), size(zMesh))
                error('AnalyticalCylinder3D:SizeMismatch', ...
                    'xMesh, yMesh, zMesh must all have identical sizes.');
            end

            [frac, insideMask] = obj.estimateImageShape(xMesh, yMesh, zMesh);

            perimMask = bwperim(insideMask);

            xOut = xMesh(perimMask);
            yOut = yMesh(perimMask);
            zOut = zMesh(perimMask);

            xOut = xOut(:);
            yOut = yOut(:);
            zOut = zOut(:);
        end
    end

    %% BODY-frame analytic FT and inside-test
    methods (Access = protected)

        % ----------- Analytic Fourier transform (BODY frame) -------------
        function S = bodyKspace(obj, kx_body, ky_body, kz_body)
            arguments
                obj
                kx_body double
                ky_body double
                kz_body double
            end

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('AnalyticalCylinder3D:SizeMismatch');
            end

            R = obj.R_mm;
            L = obj.L_mm;

            % Radial term
            k_perp = sqrt(kx_body.^2 + ky_body.^2);
            radial = zeros(size(k_perp));

            idx = (k_perp ~= 0);
            if any(idx(:))
                x_r = 2*pi*R .* k_perp(idx);
                radial(idx) = R .* besselj(1, x_r) ./ k_perp(idx);
            end
            radial(~idx) = pi * R^2;   % limit as k_perp → 0

            % Axial term
            zfac = zeros(size(kz_body));
            idx = (kz_body ~= 0);
            if any(idx(:))
                x_z = pi*L .* kz_body(idx);
                zfac(idx) = L .* (sin(x_z) ./ x_z);
            end
            zfac(~idx) = L;

            S = radial .* zfac;
        end

        % ----------------- Inside test (BODY frame) ----------------------
        function frac = percentInsideShape(obj, xb, yb, zb)
            % 0/1 mask (upgradeable later to partial volume)
            R = obj.R_mm;
            L = obj.L_mm;

            inside = (xb.^2 + yb.^2 <= R^2) & (abs(zb) <= L/2);
            frac = double(inside);
        end
    end
end
