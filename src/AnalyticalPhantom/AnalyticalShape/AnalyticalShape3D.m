classdef (Abstract) AnalyticalShape3D < handle
    % AnalyticalShape3D
    %   Base class for 3D analytic shapes in k-space.
    %
    %   Responsibilities:
    %     - Store pose (center + orientation) in WORLD coordinates.
    %     - Convert WORLD k-space coords (kx,ky,kz) [cycles/mm]
    %       into BODY-frame coords for the derived shape.
    %     - Apply translation phase ramps in WORLD coordinates.
    %     - Expose a unified `shapeChanged` event for any parameter change
    %       that affects the analytic FT (pose, size, etc.).
    %
    %   Fourier convention:
    %     S(k) = ∭ ρ(r) exp(-i 2π k·r) dV
    %
    %   BODY frame:
    %     - Defined by the derived class.
    %     - bodyKspace(kx_body, ky_body, kz_body) must implement the
    %       analytic FT of the shape in BODY coordinates, centered at 0.
    %
    %   WORLD frame:
    %     - `center` is [x0,y0,z0] in mm.
    %     - `roll_deg`, `pitch_deg`, `yaw_deg` (deg) define the body→world
    %       rotation via Rz(yaw)*Ry(pitch)*Rx(roll).

    %% Events
    events
        % Fired whenever any property that changes the shape FT is modified.
        % This includes pose changes (center/orientation) and any
        % shape-specific parameters in derived classes, when they call
        % markShapeChanged().
        shapeChanged
    end

    %% Properties
    properties
        % Center of the shape in WORLD coordinates [mm]
        center (1,3) double = [0 0 0];

        % Orientation in degrees (body → world)
        roll_deg  (1,1) double = 0;   % rotation about x-axis
        pitch_deg (1,1) double = 0;   % rotation about y-axis
        yaw_deg   (1,1) double = 0;   % rotation about z-axis
    end

    %% Constructor
    methods
        function obj = AnalyticalShape3D(center, roll_deg, pitch_deg, yaw_deg)
            % AnalyticalShape3D  Construct base shape with optional pose.
            %
            %   obj = AnalyticalShape3D()
            %   obj = AnalyticalShape3D(center)
            %   obj = AnalyticalShape3D(center, roll_deg, pitch_deg, yaw_deg)
            %
            %   center   : [1x3] double, [x0,y0,z0] in mm
            %   roll_deg : scalar double, rotation about x (deg)
            %   pitch_deg: scalar double, rotation about y (deg)
            %   yaw_deg  : scalar double, rotation about z (deg)

            if nargin >= 1 && ~isempty(center)
                obj.setCenter(center);
            end

            if nargin == 4
                obj.setOrientation(roll_deg, pitch_deg, yaw_deg);
            elseif nargin > 1 && nargin ~= 4
                error(['AnalyticalShape3D:Constructor: ', ...
                       'Provide either no orientation or all three angles.']);
            end
        end
    end

    %% Rotation helper
    methods (Access = protected)
        function R = calculateRotationMatrix(obj)
            % calculateRotationMatrix  Compute body→world rotation matrix.
            %
            %   R = calculateRotationMatrix(obj)
            %
            %   Uses intrinsic rotations:
            %     R = Rz(yaw) * Ry(pitch) * Rx(roll)

            r = deg2rad(obj.roll_deg);
            p = deg2rad(obj.pitch_deg);
            y = deg2rad(obj.yaw_deg);

            Rx = [1  0      0;
                  0  cos(r) -sin(r);
                  0  sin(r)  cos(r)];

            Ry = [ cos(p) 0 sin(p);
                   0      1 0;
                  -sin(p) 0 cos(p)];

            Rz = [cos(y) -sin(y) 0;
                  sin(y)  cos(y) 0;
                  0        0     1];

            R = Rz * Ry * Rx;   % body → world
        end

        function markShapeChanged(obj)
            % markShapeChanged  Notify listeners that shape FT has changed.
            %
            %   Call this from:
            %     - setCenter / setOrientation (already done in base class)
            %     - Any derived-class setter that changes geometry or
            %       material parameters impacting bodyKspace().
            notify(obj, 'shapeChanged');
        end
    end

    %% Pose setters
    methods
        function setCenter(obj, newCenter)
            % setCenter  Set shape center in WORLD coordinates.
            %
            %   setCenter(obj, newCenter)
            %
            %   newCenter must be a 1x3 double vector [x0,y0,z0] in mm.

            arguments
                obj
                newCenter (1,3) double
            end

            % Ensure row vector
            newCenter = double(newCenter(:)).';

            oldCenter = obj.center;
            if ~isequal(oldCenter, newCenter)
                obj.center = newCenter;
                obj.markShapeChanged();
            end
        end

        function setOrientation(obj, roll_deg, pitch_deg, yaw_deg)
            % setOrientation  Set Euler angles (body→world) in degrees.
            %
            %   setOrientation(obj, roll_deg, pitch_deg, yaw_deg)
            %
            %   roll_deg  : rotation about x-axis (deg)
            %   pitch_deg : rotation about y-axis (deg)
            %   yaw_deg   : rotation about z-axis (deg)

            arguments
                obj
                roll_deg  (1,1) double {mustBeFinite}
                pitch_deg (1,1) double {mustBeFinite}
                yaw_deg   (1,1) double {mustBeFinite}
            end

            old = [obj.roll_deg, obj.pitch_deg, obj.yaw_deg];
            new = [roll_deg,     pitch_deg,     yaw_deg    ];

            if ~isequal(old, new)
                obj.roll_deg  = roll_deg;
                obj.pitch_deg = pitch_deg;
                obj.yaw_deg   = yaw_deg;

                obj.markShapeChanged();
            end
        end
    end

    %% K-space evaluation
    methods
        function S = kspace(obj, kx, ky, kz)
            % kspace  Evaluate analytic 3D FT of the shape in WORLD k-space.
            %
            %   S = kspace(obj, kx, ky, kz)
            %
            %   Inputs:
            %     kx, ky, kz : arrays of identical size, double,
            %                  WORLD spatial frequencies [cycles/mm].
            %
            %   Output:
            %     S : complex double, same size as kx, analytic FT
            %         of the shape including rotation + translation.
            %
            %   Steps:
            %     1) WORLD → BODY k (apply inverse rotation).
            %     2) Evaluate bodyKspace(kx_body, ky_body, kz_body).
            %     3) Apply WORLD translation phase ramp exp(-i2π k·center).

            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('AnalyticalShape3D:kspace:SizeMismatch', ...
                      'kx, ky, kz must have identical size.');
            end

            % Fast path: no rotation
            if ~any([obj.roll_deg, obj.pitch_deg, obj.yaw_deg])
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                inputSize = size(kx);

                % Vectorize WORLD k
                kxv = kx(:);
                kyv = ky(:);
                kzv = kz(:);
                K_world = [kxv, kyv, kzv].';  % 3×N

                % WORLD → BODY: k_body = R^T * k_world
                R = obj.calculateRotationMatrix();
                K_body = R.' * K_world;

                % Reshape back to input size
                kxb = reshape(K_body(1,:), inputSize);
                kyb = reshape(K_body(2,:), inputSize);
                kzb = reshape(K_body(3,:), inputSize);
            end

            % 1) BODY-frame analytic FT (shape-specific)
            S_body = obj.bodyKspace(kxb, kyb, kzb);

            % 2) WORLD translation: apply phase ramp if center ≠ 0
            if any(obj.center ~= 0)
                phaseRamp = exp(-1i * 2*pi * ( ...
                    kx * obj.center(1) + ...
                    ky * obj.center(2) + ...
                    kz * obj.center(3)));
                S = S_body .* phaseRamp;
            else
                S = S_body;
            end
        end
    end

    %% Listener convenience methods
    methods
        function lh = addShapeChangedListener(obj, callback)
            % addShapeChangedListener  Attach listener for shapeChanged.
            %
            %   lh = addShapeChangedListener(obj, @(src,evt)...)
            %
            %   Fires whenever:
            %     - center changes (setCenter),
            %     - orientation changes (setOrientation),
            %     - any derived-class parameter changes and calls
            %       obj.markShapeChanged().

            arguments
                obj
                callback (1,1) function_handle
            end

            lh = addlistener(obj, 'shapeChanged', callback);
        end

        function removeListener(~, listenerHandle)
            % removeListener  Convenience wrapper to delete a listener.
            %
            %   removeListener(obj, lh)

            arguments
                ~
                listenerHandle
            end

            if ~isempty(listenerHandle) && isvalid(listenerHandle)
                delete(listenerHandle);
            end
        end
    end

    %% Abstract FT in body frame
    methods (Abstract, Access = protected)
        % bodyKspace  Analytic FT in BODY coordinates (center at origin).
        %
        %   S = bodyKspace(obj, kx_body, ky_body, kz_body)
        %
        %   Inputs:
        %     kx_body, ky_body, kz_body : BODY-frame spatial frequencies
        %                                  [cycles/mm], same size.
        %
        %   Output:
        %     S : complex double, same size as kx_body.

        S = bodyKspace(obj, kx_body, ky_body, kz_body);
    end
end
