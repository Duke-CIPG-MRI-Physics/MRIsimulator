classdef (Abstract) AnalyticalShape3D < handle & matlab.mixin.Heterogeneous
    % AnalyticalShape3D
    %   Base class for analytic 3D shapes defined in BODY coordinates,
    %   but positioned + oriented in WORLD coordinates.
    %
    %   Responsibilities:
    %     • Store pose (center + orientation)
    %     • Convert WORLD → BODY k-space coordinates
    %     • Convert WORLD → BODY spatial coordinates
    %     • Apply WORLD translation phase ramp
    %     • Store an intensity multiplier for k-space scaling
    %     • Fire a unified shapeChanged event whenever geometry changes
    %     • Provide estimateImageShape() for slice-wise shape masks
    %     • Require subclasses to implement:
    %           bodyKspace(kx_body, ky_body, kz_body)
    %           percentInsideShape(xb, yb, zb)
    %
    %   BODY frame:
    %       Shape is centered at (0,0,0), oriented with z-axis “up”.
    %
    %   WORLD frame:
    %       • center        = [x0,y0,z0] (mm)
    %       • rollPitchYaw = [roll pitch yaw] (deg)
    %
    %   Fourier convention:
    %       S(k) = ∭ ρ(r) exp(-i 2π k·r) dV
    %
    %   WORLD → BODY, k-space:
    %       k_body = Rᵀ * k_world
    %
    %   WORLD → BODY, spatial:
    %       r_body = Rᵀ * (r_world - center)

    %% Events ---------------------------------------------------------------
    events
        shapeChanged   % fire whenever *anything* affecting FT changes
    end

    %% Properties -----------------------------------------------------------
    properties (Access = private)
        center (1,3) double = [0 0 0];          % world center [mm]
        shapeIntensity double = 1;              % dimensionless multiplier (scalar or grid-sized)
        rollPitchYaw (1,3) double = [0 0 0];    % [roll pitch yaw] in deg
    end

    %% Constructor ----------------------------------------------------------
    methods
        function obj = AnalyticalShape3D(intensity, center, rollPitchYaw)
            % AnalyticalShape3D constructor.
            %
            %   obj = AnalyticalShape3D()
            %   obj = AnalyticalShape3D(intensity)
            %   obj = AnalyticalShape3D(intensity, center)
            %   obj = AnalyticalShape3D(intensity, center, rollPitchYaw)
            %
            %   intensity     : optional scaling factor (default: 1)
            %   center        : [1x3] WORLD center [mm] (default: [0 0 0])
            %   rollPitchYaw  : [1x3] Euler angles [deg] (default: [0 0 0])

            if nargin >= 1 && ~isempty(intensity)
                obj.setIntensity(intensity);
            else
                obj.setIntensity(1);
            end

            if nargin >= 2 && ~isempty(center)
                obj.setCenter(center);
            end

            if nargin >= 3 && ~isempty(rollPitchYaw)
                obj.setRollPitchYaw(rollPitchYaw);
            end
        end
    end

    %% Rotation + shared helpers -------------------------------------------
    methods (Access = protected)
        function R = calculateRotationMatrix(obj)
            % calculateRotationMatrix
            %   BODY→WORLD rotation matrix from Euler angles (deg).

            rpy = obj.getRollPitchYaw();
            r = deg2rad(rpy(1));
            p = deg2rad(rpy(2));
            y = deg2rad(rpy(3));

            Rx = [1 0 0;
                  0 cos(r) -sin(r);
                  0 sin(r)  cos(r)];

            Ry = [ cos(p) 0 sin(p);
                   0      1 0;
                  -sin(p) 0 cos(p)];

            Rz = [cos(y) -sin(y) 0;
                  sin(y)  cos(y) 0;
                  0        0     1];

            R = Rz * Ry * Rx;   % intrinsic rotations
        end

        function markShapeChanged(obj)
            % markShapeChanged
            %   Notify listeners that shape geometry/pose has changed.
            notify(obj, 'shapeChanged');
        end

        function [xb, yb, zb] = worldToBodyPoints(obj, x_world, y_world, z_world)
            % worldToBodyPoints
            %   Convert WORLD spatial coordinates (mm) to BODY coordinates.
            %
            %   [xb, yb, zb] = worldToBodyPoints(obj, x_world, y_world, z_world)
            %
            %   Inputs:
            %       x_world, y_world, z_world : same-size arrays of WORLD coords [mm]
            %
            %   Outputs:
            %       xb, yb, zb : BODY-frame coordinates [mm], same size as inputs.
            %
            %   Uses:
            %       r_body = Rᵀ * (r_world - center)

            arguments
                obj
                x_world double
                y_world double
                z_world double
            end

            if ~isequal(size(x_world), size(y_world), size(z_world))
                error('AnalyticalShape3D:worldToBodyPoints:SizeMismatch', ...
                      'x_world, y_world, z_world must have identical sizes.');
            end

            % Vectorize WORLD points: 3×N
            P_world = [x_world(:).'; y_world(:).'; z_world(:).'];

            % Subtract center (WORLD translation)
            P_shift = P_world - obj.center(:);   % 3×N

            % WORLD → BODY via Rᵀ
            R = obj.calculateRotationMatrix();
            P_body = R.' * P_shift;

            % Reshape back to original grid shape
            sz = size(x_world);
            xb = reshape(P_body(1,:), sz);
            yb = reshape(P_body(2,:), sz);
            zb = reshape(P_body(3,:), sz);
        end
    end

    %% Pose setters ---------------------------------------------------------
    methods
        function setCenter(obj, newCenter)
            % setCenter  Set WORLD center [mm].
            arguments
                obj
                newCenter (1,3) double
            end

            newCenter = double(newCenter(:)).';  % force row vector

            if ~isequal(obj.center, newCenter)
                obj.center = newCenter;
                obj.markShapeChanged();
            end
        end

        function setRollPitchYaw(obj, newRollPitchYaw)
            % setRollPitchYaw  Set BODY→WORLD Euler angles (deg) as a vector.
            arguments
                obj
                newRollPitchYaw (1,3) double {mustBeFinite}
            end

            newRollPitchYaw = double(newRollPitchYaw(:)).';

            if ~isequal(obj.rollPitchYaw, newRollPitchYaw)
                obj.rollPitchYaw = newRollPitchYaw;
                obj.markShapeChanged();
            end
        end

        function setOrientation(obj, roll_deg, pitch_deg, yaw_deg)
            % setOrientation  Set BODY→WORLD Euler angles (deg) using scalars.
            arguments
                obj
                roll_deg  (1,1) double {mustBeFinite}
                pitch_deg (1,1) double {mustBeFinite}
                yaw_deg   (1,1) double {mustBeFinite}
            end

            obj.setRollPitchYaw([roll_deg pitch_deg yaw_deg]);
        end

        function center = getCenter(obj)
            center = obj.center;
        end

        function rpy = getRollPitchYaw(obj)
            rpy = obj.rollPitchYaw;
        end
    end

    %% Intensity helpers ----------------------------------------------------
    methods
        function setIntensity(obj, newIntensity)
            arguments
                obj
                newIntensity double {mustBeFinite}
            end

            if ~isequal(obj.shapeIntensity, newIntensity)
                obj.shapeIntensity = newIntensity;
                obj.markShapeChanged();
            end
        end

        function intensity = getIntensity(obj)
            intensity = obj.shapeIntensity;
        end
    end

    %% K-space evaluation ---------------------------------------------------
    methods
        function S = kspace(obj, kx, ky, kz)
            % kspace
            %   Evaluate WORLD-frame analytic 3D FT of the shape, including
            %   the shapeIntensity scaling factor.

            intensity = obj.getIntensity();
            if ~isscalar(intensity) && ~isequal(size(intensity), size(kx))
                error('AnalyticalShape3D:kspace:IntensitySizeMismatch', ...
                    'intensity must be scalar or match size of kx/ky/kz.');
            end

            S = intensity .* obj.kspace_shapeOnly(kx, ky, kz);
        end
    end

    methods (Static, Access = protected)
        function paramOut = requireScalarOrSize(param, template, paramName)
            % requireScalarOrSize
            %   Enforce that a shape parameter is either scalar or matches the
            %   size of a template array (typically kx/ky/kz or xb/yb/zb).

            if isscalar(param)
                % Expand scalar parameters to match the template size so that
                % subsequent logical indexing (e.g., param(idx)) never
                % produces an out-of-bounds error when the template is larger
                % than a single element.
                paramOut = param .* ones(size(template));
            elseif isequal(size(param), size(template))
                paramOut = param;
            else
                error('AnalyticalShape3D:ParameterSizeMismatch', ...
                    'Parameter %s must be scalar or match template size.', paramName);
            end
        end
    end

    methods
        function S = kspace_shapeOnly(obj, kx, ky, kz)
            % kspace_shapeOnly
            %   Evaluate WORLD-frame analytic 3D FT of the shape geometry
            %   (no intensity scaling).
            %
            %   S = kspace_shapeOnly(obj, kx, ky, kz)
            %
            %   Inputs:
            %       kx,ky,kz : WORLD freqs [cycles/mm], same size.
            %
            %   Output:
            %       S : complex double, same size as kx.
            %
            %   Steps:
            %       1) WORLD → BODY k
            %       2) Evaluate bodyKspace() in BODY frame
            %       3) Apply WORLD translation phase ramp

            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('AnalyticalShape3D:kspace:SizeMismatch', ...
                      'kx, ky, kz must have identical sizes.');
            end

            rpy = obj.getRollPitchYaw();
            noRotation = ~any(rpy);

            if noRotation
                % Fast path: BODY and WORLD frames align
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                % WORLD → BODY transform
                inputSize = size(kx);
                K_world = [kx(:) ky(:) kz(:)].';
                R = obj.calculateRotationMatrix();
                K_body = R.' * K_world;

                kxb = reshape(K_body(1,:), inputSize);
                kyb = reshape(K_body(2,:), inputSize);
                kzb = reshape(K_body(3,:), inputSize);
            end

            % Analytic FT in BODY frame
            S_body = obj.bodyKspace(kxb, kyb, kzb);

            % WORLD translation phase
            if any(obj.center ~= 0)
                c = obj.getCenter();
                phase = exp(-1i * 2*pi * ( ...
                        kx * c(1) + ...
                        ky * c(2) + ...
                        kz * c(3)));
                S = S_body .* phase;
            else
                S = S_body;
            end
        end
    end

    %% Image-domain shape estimation ---------------------------------------
    methods
        function imageShape = estimateImageShape(obj, xMesh, yMesh, zMesh)
            % estimateImageShape
            %   Compute a binary mask of the shape on an arbitrary slice plane.
            %
            %   imageShape = estimateImageShape(obj, xMesh, yMesh, zMesh)
            %
            %   Inputs:
            %       xMesh, yMesh, zMesh : 2D arrays (same size), WORLD
            %                             coordinates [mm] of pixel centers
            %                             on some slice plane (axial, coronal,
            %                             sagittal, or oblique).
            %
            %   Output:
            %       imageShape : logical mask, same size as xMesh
            %                    true where the shape intersects the pixel
            %                    (currently: "inside at pixel center").
            %
            %   Notes:
            %       • This is currently a 0/1 mask based on the BODY-frame
            %         point-in-shape test at each pixel center.
            %       • Later, percentInsideShape() can be extended to return
            %         partial-volume fractions if desired.

            arguments
                obj
                xMesh double
                yMesh double
                zMesh double
            end

            if ~isequal(size(xMesh), size(yMesh), size(zMesh))
                error('AnalyticalShape3D:estimateImageShape:SizeMismatch', ...
                      'xMesh, yMesh, zMesh must have identical sizes.');
            end

            % WORLD → BODY coordinates for all pixels
            [xb, yb, zb] = obj.worldToBodyPoints(xMesh, yMesh, zMesh);

            % Subclass computes inside/outside mask in BODY frame
            mask = obj.percentInsideShape(xb, yb, zb);

            % Basic checks
            if ~isequal(size(mask), size(xMesh))
                error('AnalyticalShape3D:estimateImageShape:SizeMismatch', ...
                      'percentInsideShape must return same size as xb,yb,zb.');
            end

            % Coerce to logical 0/1 mask
            imageShape = (mask ~= 0);
        end

        function image = estimateImage(obj, xMesh, yMesh, zMesh)
            % estimateImage
            %   Estimate image intensity by masking the shape and scaling by
            %   the shape's intensity.

            image = obj.getIntensity() .* obj.estimateImageShape(xMesh, yMesh, zMesh);
        end
    end

    %% Listener helpers -----------------------------------------------------
    methods
        function lh = addShapeChangedListener(obj, callback)
            % addShapeChangedListener
            %   Attach listener to the shapeChanged event.
            %
            %   lh = addShapeChangedListener(obj, @(src,evt)...)
            lh = addlistener(obj, 'shapeChanged', callback);
        end

        function removeListener(~, lh)
            % removeListener
            %   Convenience wrapper to delete a listener handle.
            if ~isempty(lh) && isvalid(lh)
                delete(lh);
            end
        end
    end

    %% Abstract analytic FT + inside test ----------------------------------
    methods (Abstract, Access = protected)
        % bodyKspace
        %   Analytic FT of the shape in the BODY frame.
        %
        %   S = bodyKspace(obj, kx_body, ky_body, kz_body)
        %
        %   Inputs:
        %       kx_body, ky_body, kz_body : BODY freqs [cycles/mm], same size
        %
        %   Output:
        %       S : complex double, same size as kx_body
        S = bodyKspace(obj, kx_body, ky_body, kz_body);

        % percentInsideShape
        %   BODY-frame inside/outside test on a spatial grid.
        %
        %   imageShape = percentInsideShape(obj, xb, yb, zb)
        %
        %   Inputs:
        %       xb,yb,zb : BODY coordinates [mm] of pixel centers, same size
        %
        %   Output:
        %       imageShape : numeric or logical array, same size as xb
        %                    currently expected to encode 0 (outside) and
        %                    nonzero (inside). Base class converts to logical.
        %
        %   Notes:
        %       • For now, subclasses typically implement this as a 0/1
        %         mask based on center-point inclusion tests.
        %       • In the future, this can be generalized to partial-volume
        %         fractions in [0,1], but estimateImageShape currently
        %         binarizes to inside/outside.
        imageShape = percentInsideShape(obj, xb, yb, zb);
    end

end
