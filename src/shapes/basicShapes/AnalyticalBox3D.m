classdef AnalyticalBox3D < AnalyticalShape3D
    % AnalyticalBox3D
    %   Rectangular prism aligned with BODY axes.

    properties (Access = protected)
        Lx_mm double {mustBeNonnegative} = 1;
        Ly_mm double {mustBeNonnegative} = 1;
        Lz_mm double {mustBeNonnegative} = 1;
    end

    methods
        function obj = AnalyticalBox3D(Lx_mm, Ly_mm, Lz_mm, intensity, center, rollPitchYaw)
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
            if nargin >= 1 && ~isempty(Lx_mm)
                obj.Lx_mm = Lx_mm;
            end
            if nargin >= 2 && ~isempty(Ly_mm)
                obj.Ly_mm = Ly_mm;
            end
            if nargin >= 3 && ~isempty(Lz_mm)
                obj.Lz_mm = Lz_mm;
            end
        end

        function setSideLengths(obj, Lx_mm, Ly_mm, Lz_mm)
            obj.Lx_mm = Lx_mm;
            obj.Ly_mm = Ly_mm;
            obj.Lz_mm = Lz_mm;
        end

        function vol = calculateVolume(obj)
            vol = obj.Lx_mm .* obj.Ly_mm .* obj.Lz_mm;
        end
    end

    methods (Access = protected)
        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            S_body = obj.calculateVolume() .* ...
                sinc(kx .* obj.Lx_mm) .* ...
                sinc(ky .* obj.Ly_mm) .* ...
                sinc(kz .* obj.Lz_mm);
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            inX = abs(xb) <= obj.Lx_mm ./ 2;
            inY = abs(yb) <= obj.Ly_mm ./ 2;
            inZ = abs(zb) <= obj.Lz_mm ./ 2;
            percent = double(inX & inY & inZ);
        end
    end
end
