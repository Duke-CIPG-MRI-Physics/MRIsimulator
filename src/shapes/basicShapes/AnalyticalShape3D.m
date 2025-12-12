classdef (Abstract) AnalyticalShape3D < handle & matlab.mixin.Heterogeneous
    % AnalyticalShape3D
    %   Abstract base class for analytic 3D shapes with configurable
    %   intensity, pose, and time-varying geometry. Shape parameters can be
    %   set with static structs or a function handle that produces a
    %   parameter struct when evaluated (e.g., with time). Subclasses
    %   implement the BODY-frame Fourier transform (kspaceBaseShape) and
    %   voxel occupancy (percentInsideBody) for their specific geometries
    %   while mapping their dimensions into the common setShapeParameters
    %   /getShapeParameters API.

    properties (Access = protected)
        shapeIntensity (1,1) double = 1;
        shapeParameters = struct();
    end

    methods
        function obj = AnalyticalShape3D(intensity, shapeParameters)
            if nargin < 1 || isempty(intensity)
                intensity = 1;
            end
            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = AnalyticalShape3D.ensurePoseFields(struct());
            end

            obj.shapeIntensity = intensity;
            obj.setShapeParameters(shapeParameters);
        end

        function setShapeParameters(obj, shapeParameters)
            % setShapeParameters  Store shape parameters or a generating function handle.
            %   setShapeParameters(obj, shapeParameters) validates the provided
            %   struct or the struct returned by a function handle before
            %   storing the original input. Validation is split into required
            %   field checks (validateParameterFields) and size compatibility
            %   checks (validateParameterSize).
            paramsForValidation = obj.evaluateShapeParameters(shapeParameters, {});
            obj.validateParameterFields(paramsForValidation);
            obj.validateParameterSize(paramsForValidation);

            obj.shapeParameters = shapeParameters;
        end

        function params = getShapeParameters(obj, varargin)
            % getShapeParameters  Retrieve the current shape parameters.
            %   params = getShapeParameters(obj, ...) returns the stored struct
            %   directly or evaluates the stored function handle with any
            %   provided arguments before validation.
            params = obj.evaluateShapeParameters(obj.shapeParameters, varargin);

            obj.validateParameterFields(params);
            obj.validateParameterSize(params, varargin{:});
        end

        function setIntensity(obj, newIntensity)
            obj.shapeIntensity = newIntensity;
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
            %   Applies rotation from roll/pitch/yaw, translation to obj.pose,
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

            params = obj.getShapeParameters();
            pose = obj.extractPose(params, kx);
            rpy = [pose.roll_deg, pose.pitch_deg, pose.yaw_deg];
            noRotation = ~any(rpy(:));

            if noRotation
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                [kxb, kyb, kzb] = obj.rotateWorldToBody(kx, ky, kz, pose.roll_deg, pose.pitch_deg, pose.yaw_deg);
            end

            S_body = obj.kspaceBaseShape(kxb, kyb, kzb);

            if any(pose.centerX(:) ~= 0) || any(pose.centerY(:) ~= 0) || any(pose.centerZ(:) ~= 0)
                phase = exp(-1i * 2*pi * ( ...
                        kx .* pose.centerX + ...
                        ky .* pose.centerY + ...
                        kz .* pose.centerZ));
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

            params = obj.getShapeParameters();
            pose = obj.extractPose(params, x);

            xBody = x - pose.centerX;
            yBody = y - pose.centerY;
            zBody = z - pose.centerZ;

            if any([pose.roll_deg(:); pose.pitch_deg(:); pose.yaw_deg(:)])
                [xBody, yBody, zBody] = obj.rotateWorldToBody(xBody, yBody, zBody, ...
                    pose.roll_deg, pose.pitch_deg, pose.yaw_deg);
            end

            percent = obj.percentInsideBody(xBody, yBody, zBody);
        end
    end

    methods (Access = protected)
        function params = evaluateShapeParameters(~, parameterSource, args)
            if isa(parameterSource, 'function_handle')
                params = parameterSource(args{:});
            else
                params = parameterSource;
            end

            if ~isstruct(params)
                error('AnalyticalShape3D:ShapeParameters:InvalidType', ...
                    'Shape parameters must be provided as a struct.');
            end
        end

        function pose = extractPose(obj, params, referenceArray, varargin)
            poseX = obj.evaluateParameterValue(params.pose.center.x_mm, varargin{:});
            poseY = obj.evaluateParameterValue(params.pose.center.y_mm, varargin{:});
            poseZ = obj.evaluateParameterValue(params.pose.center.z_mm, varargin{:});
            roll = obj.evaluateParameterValue(params.pose.roll_deg, varargin{:});
            pitch = obj.evaluateParameterValue(params.pose.pitch_deg, varargin{:});
            yaw = obj.evaluateParameterValue(params.pose.yaw_deg, varargin{:});

            pose.centerX = obj.requireScalarOrSize(poseX, referenceArray, 'poseCenterX');
            pose.centerY = obj.requireScalarOrSize(poseY, referenceArray, 'poseCenterY');
            pose.centerZ = obj.requireScalarOrSize(poseZ, referenceArray, 'poseCenterZ');
            pose.roll_deg = obj.requireScalarOrSize(roll, referenceArray, 'poseRoll');
            pose.pitch_deg = obj.requireScalarOrSize(pitch, referenceArray, 'posePitch');
            pose.yaw_deg = obj.requireScalarOrSize(yaw, referenceArray, 'poseYaw');
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

        function validateParameterFields(obj, params)
            if ~isfield(params, 'pose') || ~isstruct(params.pose)
                error('AnalyticalShape3D:ShapeParameters:MissingPose', ...
                    'shapeParameters.pose must be provided as a struct.');
            end

            requiredCenter = {'x_mm', 'y_mm', 'z_mm'};
            if ~isfield(params.pose, 'center') || ~isstruct(params.pose.center)
                error('AnalyticalShape3D:ShapeParameters:MissingCenter', ...
                    'shapeParameters.pose.center must be provided as a struct.');
            end

            for idx = 1:numel(requiredCenter)
                if ~isfield(params.pose.center, requiredCenter{idx})
                    error('AnalyticalShape3D:ShapeParameters:MissingCenterField', ...
                        'shapeParameters.pose.center.%s is required.', requiredCenter{idx});
                end
            end

            requiredPose = {'roll_deg', 'pitch_deg', 'yaw_deg'};
            for idx = 1:numel(requiredPose)
                if ~isfield(params.pose, requiredPose{idx})
                    error('AnalyticalShape3D:ShapeParameters:MissingPoseField', ...
                        'shapeParameters.pose.%s is required.', requiredPose{idx});
                end
            end
        end

        function validateParameterSize(obj, params, varargin)
            values = obj.collectParameterValues(params, varargin{:});

            matrixSize = [];
            for idx = 1:numel(values)
                value = values{idx};
                numericValue = obj.evaluateParameterValue(value, varargin{:});
                validateattributes(numericValue, {'numeric'}, {'nonempty', 'real', 'finite'});

                if isscalar(numericValue)
                    continue;
                end

                currentSize = size(numericValue);
                if isempty(matrixSize)
                    matrixSize = currentSize;
                elseif ~isequal(matrixSize, currentSize)
                    error('AnalyticalShape3D:ShapeParameters:SizeMismatch', ...
                        'All non-scalar parameters must share the same size.');
                end
            end
        end

        function values = collectParameterValues(obj, params, varargin)
            values = {};
            fields = fieldnames(params);
            for idx = 1:numel(fields)
                fieldValue = params.(fields{idx});
                if isstruct(fieldValue)
                    nestedValues = obj.collectParameterValues(fieldValue, varargin{:});
                    values = [values, nestedValues]; %#ok<AGROW>
                else
                    values{end+1} = fieldValue; %#ok<AGROW>
                end
            end
        end

        function value = evaluateParameterValue(~, valueOrHandle, varargin)
            if isa(valueOrHandle, 'function_handle')
                value = valueOrHandle(varargin{:});
            else
                value = valueOrHandle;
            end
        end

        function [xb, yb, zb] = rotateWorldToBody(~, x, y, z, roll_deg, pitch_deg, yaw_deg)
            roll_rad = deg2rad(roll_deg);
            pitch_rad = deg2rad(pitch_deg);
            yaw_rad = deg2rad(yaw_deg);

            cr = cos(roll_rad); sr = sin(roll_rad);
            cp = cos(pitch_rad); sp = sin(pitch_rad);
            cy = cos(yaw_rad); sy = sin(yaw_rad);

            xb = cy .* cp .* x + sy .* cp .* y - sp .* z;
            yb = (cy .* sp .* sr - sy .* cr) .* x + (sy .* sp .* sr + cy .* cr) .* y + cp .* sr .* z;
            zb = (cy .* sp .* cr + sy .* sr) .* x + (sy .* sp .* cr - cy .* sr) .* y + cp .* cr .* z;
        end
    end

    methods (Static, Access = protected)
        function params = ensurePoseFields(params)
            poseDefaults = struct( ...
                'center', struct( ...
                    'x_mm', 0, ...
                    'y_mm', 0, ...
                    'z_mm', 0), ...
                'roll_deg', 0, ...
                'pitch_deg', 0, ...
                'yaw_deg', 0);

            if ~isfield(params, 'pose') || ~isstruct(params.pose)
                params.pose = poseDefaults;
                return;
            end

            if ~isfield(params.pose, 'center') || ~isstruct(params.pose.center)
                params.pose.center = poseDefaults.center;
            else
                if ~isfield(params.pose.center, 'x_mm')
                    params.pose.center.x_mm = poseDefaults.center.x_mm;
                end
                if ~isfield(params.pose.center, 'y_mm')
                    params.pose.center.y_mm = poseDefaults.center.y_mm;
                end
                if ~isfield(params.pose.center, 'z_mm')
                    params.pose.center.z_mm = poseDefaults.center.z_mm;
                end
            end

            if ~isfield(params.pose, 'roll_deg')
                params.pose.roll_deg = poseDefaults.roll_deg;
            end
            if ~isfield(params.pose, 'pitch_deg')
                params.pose.pitch_deg = poseDefaults.pitch_deg;
            end
            if ~isfield(params.pose, 'yaw_deg')
                params.pose.yaw_deg = poseDefaults.yaw_deg;
            end
        end
    end

    methods (Abstract, Access = protected)
        % kspaceBaseShape  BODY-frame k-space (no rotation/translation/intensity).
        S_body = kspaceBaseShape(obj, kx_body, ky_body, kz_body)
        percent = percentInsideBody(obj, xb, yb, zb)
    end
end
