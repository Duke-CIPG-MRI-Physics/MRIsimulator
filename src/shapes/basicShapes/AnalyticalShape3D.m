classdef (Abstract) AnalyticalShape3D < handle & matlab.mixin.Heterogeneous
    % AnalyticalShape3D
    %   Abstract base class for analytic 3D shapes with configurable
    %   intensity, translation, orientation, and time-varying geometry. Shape
    %   parameters can be set with static scalars/vectors or a function handle
    %   that produces a parameter struct when evaluated (e.g., with time).
    %   Subclasses implement the BODY-frame Fourier transform (kspaceBaseShape)
    %   and voxel occupancy (percentInsideBody) for their specific geometries
    %   while mapping their dimensions into the common setShapeParameters
    %   /getShapeParameters API.
    %
    %   Fourier transform workflow:
    %     kspaceBaseShape        – BODY-frame, unshifted/unrotated, no intensity.
    %     kspaceWorldPlacedShape – WORLD-frame, rotated/translated, no intensity.
    %     kspace                 – WORLD-frame, rotated/translated with intensity
    %                               applied to the centered and oriented shape.

    properties (Access = protected)
        shapeIntensity (1,1) double = 1;
        center (1,3) double = [0 0 0];
        rollPitchYaw_deg (1,3) double = [0 0 0];
        shapeParameters = struct();
    end

    methods
        function obj = AnalyticalShape3D(intensity, center, rollPitchYaw)
            if nargin >= 1 && ~isempty(intensity)
                obj.shapeIntensity = intensity;
            end
            if nargin >= 2 && ~isempty(center)
                obj.center = center;
            end
            if nargin >= 3 && ~isempty(rollPitchYaw)
                obj.rollPitchYaw_deg = rollPitchYaw;
            end
        end


        function setShapeParameters(obj, shapeParameters)
            % setShapeParameters  Store shape parameters or a generating function handle.
            %   setShapeParameters(obj, shapeParameters) validates a struct
            %   definition via validateParameters and stores it, or stores a
            %   function handle that returns such a struct when invoked (e.g.,
            %   time-varying geometries).
            if isa(shapeParameters, 'function_handle')
                obj.shapeParameters = shapeParameters;
                return;
            end

            % Validate that the parameters are valid
            obj.shapeParameters = obj.validateParameters(shapeParameters);
        end

        function params = getShapeParameters(obj, varargin)
            % getShapeParameters  Retrieve the current shape parameters.
            %   params = getShapeParameters(obj, ...) returns the stored struct
            %   directly or evaluates the stored function handle with any
            %   provided arguments before validation.
            if isa(obj.shapeParameters, 'function_handle')
                params = obj.shapeParameters(varargin{:});
            else
                params = obj.shapeParameters;
            end

            params = obj.validateParameters(params);
        end

        function setIntensity(obj, newIntensity)
            obj.shapeIntensity = newIntensity;
        end

        function setCenter(obj, newCenter)
            obj.center = newCenter;
        end

        function setOrientation(obj, roll_deg, pitch_deg, yaw_deg)
            obj.rollPitchYaw_deg = [roll_deg, pitch_deg, yaw_deg];
        end

        function setRollPitchYaw(obj, rollPitchYaw)
            obj.rollPitchYaw_deg = rollPitchYaw;
        end

        function c = getCenter(obj)
            c = obj.center;
        end

        function rpy = getRollPitchYaw(obj)
            rpy = obj.rollPitchYaw_deg;
        end

        function vol = estimateImage(obj, xMesh, yMesh, zMesh)
            arguments
                obj
                xMesh double
                yMesh double
                zMesh double
            end

            percent = obj.percentInsideShape(xMesh, yMesh, zMesh);
            vol = percent .* obj.shapeIntensity;
        end

        function S = kspace(obj, kx, ky, kz)
            % kspace  WORLD-frame k-space including center, orientation, and intensity.
            %   S = kspace(obj, kx, ky, kz)
            %   Applies rotation from roll/pitch/yaw, translation to obj.center,
            %   and then intensity scaling for the posed shape.
            arguments
                obj
                kx double
                ky double
                kz double
            end

            S_shape = obj.kspaceWorldPlacedShape(kx, ky, kz);
            S = S_shape .* obj.shapeIntensity;
        end

        function S = kspaceWorldPlacedShape(obj, kx, ky, kz)
            % kspaceWorldPlacedShape  WORLD-frame k-space of posed shape (no intensity).
            %   S = kspaceWorldPlacedShape(obj, kx, ky, kz)
            %   Applies rotation and translation to the BODY-frame k-space.
            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('AnalyticalShape3D:kspaceWorldPlacedShape:SizeMismatch', ...
                    'kx, ky, kz must have identical sizes.');
            end

            rpy = obj.getRollPitchYaw();
            noRotation = ~any(rpy);

            if noRotation
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                inputSize = size(kx);
                K_world = [kx(:) ky(:) kz(:)].';
                R = obj.calculateRotationMatrix();
                K_body = R.' * K_world;

                kxb = reshape(K_body(1,:), inputSize);
                kyb = reshape(K_body(2,:), inputSize);
                kzb = reshape(K_body(3,:), inputSize);
            end

            S_body = obj.kspaceBaseShape(kxb, kyb, kzb);

            c = obj.getCenter();
            if any(c(:) ~= 0)
                if size(c, 2) ~= 3
                    error('AnalyticalShape3D:kspaceWorldPlacedShape:CenterSizeMismatch', ...
                        'Center must have 3 columns for x, y, z.');
                end

                cx = obj.requireScalarOrSize(c(:,1), kx, 'centerX');
                cy = obj.requireScalarOrSize(c(:,2), ky, 'centerY');
                cz = obj.requireScalarOrSize(c(:,3), kz, 'centerZ');

                phase = exp(-1i * 2*pi * ( ...
                        kx .* cx + ...
                        ky .* cy + ...
                        kz .* cz));
                S = S_body .* phase;
            else
                S = S_body;
            end
        end

        function percent = percentInsideShape(obj, x, y, z)
            arguments
                obj
                x double
                y double
                z double
            end

            if ~isequal(size(x), size(y), size(z))
                error('AnalyticalShape3D:percentInsideShape:SizeMismatch', ...
                    'x, y, z must have identical sizes.');
            end

            c = obj.getCenter();
            rpy = obj.getRollPitchYaw();

            xBody = x - c(1);
            yBody = y - c(2);
            zBody = z - c(3);

            if any(rpy)
                coords = [xBody(:) yBody(:) zBody(:)].';
                R = obj.calculateRotationMatrix();
                bodyCoords = R.' * coords;
                inputSize = size(xBody);
                xBody = reshape(bodyCoords(1,:), inputSize);
                yBody = reshape(bodyCoords(2,:), inputSize);
                zBody = reshape(bodyCoords(3,:), inputSize);
            end

            percent = obj.percentInsideBody(xBody, yBody, zBody);
        end
    end

    methods (Access = protected)
        function R = calculateRotationMatrix(obj)
            r = deg2rad(obj.rollPitchYaw_deg(1));
            p = deg2rad(obj.rollPitchYaw_deg(2));
            y = deg2rad(obj.rollPitchYaw_deg(3));

            Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
            Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
            Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];

            R = Rz * Ry * Rx;
        end

        function vec = requireScalarOrSize(~, value, reference, id)
            if isscalar(value)
                vec = value;
            elseif isequal(size(value), size(reference))
                vec = value;
            else
                error(sprintf('AnalyticalShape3D:%s:SizeMismatch', id), ...
                    'Values must be scalar or match the reference array size.');
            end
        end
        
    end

    methods (Abstract, Access = protected)
        % kspaceBaseShape  BODY-frame k-space (no rotation/translation/intensity).
        S_body = kspaceBaseShape(obj, kx_body, ky_body, kz_body)
        percent = percentInsideBody(obj, xb, yb, zb)
        params = validateParameters(obj, params)
    end
end
