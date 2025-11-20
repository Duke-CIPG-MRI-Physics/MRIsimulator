classdef AnalyticalBox3D < AnalyticalShape3D
    % AnalyticalBox3D
    %   Solid axis-aligned rectangular box with analytic 3D FT.
    %
    %   Body-frame definition:
    %     - x ∈ [-Lx/2, +Lx/2]
    %     - y ∈ [-Ly/2, +Ly/2]
    %     - z ∈ [-Lz/2, +Lz/2]
    %
    %   World-frame pose handled by AnalyticalShape3D:
    %     - center: [x0,y0,z0] in mm
    %     - roll_deg, pitch_deg, yaw_deg: orientation (deg), body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)
    
    properties
        Lx (1,1) double = 10;   % box length in x (mm)
        Ly (1,1) double = 10;   % box length in y (mm)
        Lz (1,1) double = 10;   % box length in z (mm)
    end
    
    methods
        function obj = AnalyticalBox3D(center, roll_deg, pitch_deg, yaw_deg, ...
                                       Lx, Ly, Lz)
            % Constructor
            %
            %   obj = AnalyticalBox3D(center, roll, pitch, yaw, Lx, Ly, Lz)
            %
            %   center: [x0,y0,z0] in mm (1x3)
            %   roll/pitch/yaw: orientation in deg (body -> world)
            %   Lx, Ly, Lz: side lengths in mm
            
            obj@AnalyticalShape3D(center, roll_deg, pitch_deg, yaw_deg);
            
            if nargin >= 5 && ~isempty(Lx), obj.Lx = Lx; end
            if nargin >= 6 && ~isempty(Ly), obj.Ly = Ly; end
            if nargin >= 7 && ~isempty(Lz), obj.Lz = Lz; end
        end
    end
    
    methods (Access = protected)
        function S = bodyKspace(obj, kx_body, ky_body, kz_body)
            % bodyKspace
            %   Analytic FT of a solid box centered at origin in body frame.
            %
            %   Uses MATLAB's sinc(u) = sin(pi*u)/(pi*u).
            %   For 1D box of length L:
            %       F(k) = L * sinc(k*L)
            %
            %   So in 3D:
            %       F(kx,ky,kz) = Lx*Ly*Lz *
            %                     sinc(kx*Lx) * sinc(ky*Ly) * sinc(kz*Lz)
            
            Lx = obj.Lx;
            Ly = obj.Ly;
            Lz = obj.Lz;
            
            Sx = Lx .* sinc(kx_body .* Lx);
            Sy = Ly .* sinc(ky_body .* Ly);
            Sz = Lz .* sinc(kz_body .* Lz);
            
            S = Sx .* Sy .* Sz;
        end
    end
end
