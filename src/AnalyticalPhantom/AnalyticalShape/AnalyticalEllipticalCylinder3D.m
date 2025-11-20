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
    %     - roll_deg, pitch_deg, yaw_deg: orientation (deg), body -> world
    %
    %   kx, ky, kz in cycles/mm; phase term exp(-i 2*pi * k · r)
    
    properties
        a_mm (1,1) double = 5;    % semi-axis in x (mm)
        b_mm (1,1) double = 3;    % semi-axis in y (mm)
        L_mm (1,1) double = 20;   % length in z (mm)
    end
    
    methods
        function obj = AnalyticalEllipticalCylinder3D(center, roll_deg, pitch_deg, yaw_deg, ...
                                                     a_mm, b_mm, L_mm)
            % Constructor
            %
            %   obj = AnalyticalEllipticalCylinder3D(center, roll, pitch, yaw, ...
            %                                        a_mm, b_mm, L_mm)
            %
            %   center: [x0,y0,z0] in mm (1x3)
            %   roll/pitch/yaw: orientation in deg (body -> world)
            %   a_mm, b_mm: semi-axes of ellipse in x,y (mm)
            %   L_mm: length in mm
            
            obj@AnalyticalShape3D(center, roll_deg, pitch_deg, yaw_deg);
            
            if nargin >= 5 && ~isempty(a_mm), obj.a_mm = a_mm; end
            if nargin >= 6 && ~isempty(b_mm), obj.b_mm = b_mm; end
            if nargin >= 7 && ~isempty(L_mm), obj.L_mm = L_mm; end
        end
    end
    
    methods (Access = protected)
        function S = bodyKspace(obj, kx_body, ky_body, kz_body)
            % bodyKspace
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
            
            a = obj.a_mm;
            b = obj.b_mm;
            L = obj.L_mm;
            
            % Scaled radial frequency
            k_perp_prime = sqrt( (a * kx_body).^2 + (b * ky_body).^2 );
            
            % Radial part (base cylinder radius R0 = 1, then scale by a*b)
            radial = zeros(size(k_perp_prime));
            idx_r = (k_perp_prime ~= 0);
            if any(idx_r(:))
                x_r = 2*pi * k_perp_prime(idx_r);    % R0=1
                radial(idx_r) = besselj(1, x_r) ./ k_perp_prime(idx_r);
            end
            % k_perp' -> 0 limit: area of ellipse = pi*a*b
            radial(~idx_r) = pi * a * b;
            
            % Axial part (same as finite cylinder of length L)
            kz = kz_body;
            zfac = zeros(size(kz));
            idx_z = (kz ~= 0);
            if any(idx_z(:))
                x_z = pi * L .* kz(idx_z);
                zfac(idx_z) = L .* (sin(x_z) ./ x_z);   % L * sinc(pi*L*kz)
            end
            % kz -> 0 limit: length
            zfac(~idx_z) = L;
            
            S = radial .* zfac;
        end
    end
end
