classdef (Abstract) AnalyticalShape3D < handle
    % AnalyticalShape3D
    %   Abstract base class for analytic 3D shapes with configurable
    %   intensity, translation, and orientation. Subclasses implement the
    %   BODY-frame Fourier transform (bodyKspace) and voxel occupancy
    %   (percentInsideBody) for their specific geometries.

    events
        shapeChanged
    end

    properties (Access = protected)
        shapeIntensity (1,1) double = 1;
        center (1,3) double = [0 0 0];
        rollPitchYaw_deg (1,3) double = [0 0 0];
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

        function setIntensity(obj, newIntensity)
            obj.shapeIntensity = newIntensity;
            obj.markShapeChanged();
        end

        function setCenter(obj, newCenter)
            obj.center = newCenter;
            obj.markShapeChanged();
        end

        function setOrientation(obj, roll_deg, pitch_deg, yaw_deg)
            obj.rollPitchYaw_deg = [roll_deg, pitch_deg, yaw_deg];
            obj.markShapeChanged();
        end

        function setRollPitchYaw(obj, rollPitchYaw)
            obj.rollPitchYaw_deg = rollPitchYaw;
            obj.markShapeChanged();
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
            arguments
                obj
                kx double
                ky double
                kz double
            end

            S_body = obj.kspace_shapeOnly(kx, ky, kz);
            S = S_body .* obj.shapeIntensity;
        end

        function S = kspace_shapeOnly(obj, kx, ky, kz)
            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('AnalyticalShape3D:kspace_shapeOnly:SizeMismatch', ...
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

            S_body = obj.bodyKspace(kxb, kyb, kzb);

            c = obj.getCenter();
            if any(c(:) ~= 0)
                if size(c, 2) ~= 3
                    error('AnalyticalShape3D:kspace_shapeOnly:CenterSizeMismatch', ...
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

        function listenerHandle = addShapeChangedListener(obj, callback)
            listenerHandle = addlistener(obj, 'shapeChanged', callback);
        end
    end

    methods (Access = protected)
        function markShapeChanged(obj)
            notify(obj, 'shapeChanged');
        end

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

    methods (Abstract)
        vol = calculateVolume(obj)
    end

    methods (Abstract, Access = protected)
        S_body = bodyKspace(obj, kx_body, ky_body, kz_body)
        percent = percentInsideBody(obj, xb, yb, zb)
    end
end
