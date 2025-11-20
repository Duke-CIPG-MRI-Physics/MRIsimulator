classdef AnalyticalEllipsoid3D < AnalyticalShape3D
    % AnalyticalEllipsoid3D
    %   Solid ellipsoid with analytic 3D FT in k-space.
    %
    %   Body-frame definition:
    %     (x/a_mm)^2 + (y/b_mm)^2 + (z/c_mm)^2 <= 1
    %
    %   World-frame pose handled by AnalyticalShape3D:
    %     - center: [x0,y0,z0] in mm
    %     - roll_deg, pitch_deg, yaw_deg: orientation (deg), body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)
    
    properties
        a_mm (1,1) double = 5;   % semi-axis in x (mm)
        b_mm (1,1) double = 5;   % semi-axis in y (mm)
        c_mm (1,1) double = 5;   % semi-axis in z (mm)
    end
    
    methods
        function obj = AnalyticalEllipsoid3D(center, roll_deg, pitch_deg, yaw_deg, ...
                                             a_mm, b_mm, c_mm)
            % Constructor
            %
            %   obj = AnalyticalEllipsoid3D(center, roll, pitch, yaw, ...
            %                               a_mm, b_mm, c_mm)
            %
            %   center: [x0,y0,z0] in mm (1x3)
            %   roll/pitch/yaw: orientation in deg (body -> world)
            %   a_mm, b_mm, c_mm: semi-axes in mm
            
            obj@AnalyticalShape3D(center, roll_deg, pitch_deg, yaw_deg);
            
            if nargin >= 5 && ~isempty(a_mm), obj.a_mm = a_mm; end
            if nargin >= 6 && ~isempty(b_mm), obj.b_mm = b_mm; end
            if nargin >= 7 && ~isempty(c_mm), obj.c_mm = c_mm; end
        end
    end
    
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
    end
end
