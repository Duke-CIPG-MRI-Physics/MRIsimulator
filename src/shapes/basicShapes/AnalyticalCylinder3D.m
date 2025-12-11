classdef AnalyticalCylinder3D < AnalyticalShape3D
    % AnalyticalCylinder3D
    %   Finite cylinder aligned with the BODY +z axis. Radius/length can be
    %   static scalars/vectors, a struct of those values, or supplied by a
    %   function handle returning such a struct for time-varying geometries.

    methods
        function obj = AnalyticalCylinder3D(shapeParameters, intensity, center, rollPitchYaw)
            if nargin < 4 || isempty(rollPitchYaw)
                rollPitchYaw = [0 0 0];
            end
            if nargin < 3 || isempty(center)
                center = [0 0 0];
            end
            if nargin < 2 || isempty(intensity)
                intensity = 1;
            end
             
            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);

            if nargin < 1 || isempty(shapeParameters)
                shapeParameters = obj.defaultCylinderParameters();
            end
            obj.setShapeParameters(shapeParameters);
        end
    end

    methods (Access = protected)
        function params = validateParameters(~, params)
            if isa(params, 'function_handle')
                params = params();
            end

            if ~isstruct(params)
                error('AnalyticalCylinder3D:ShapeParameters:InvalidType', ...
                    'Shape parameters must be provided as a struct.');
            end

            required = {'radius_mm', 'length_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalCylinder3D:ShapeParameters:MissingField', ...
                        'Field %s is required for cylinder dimensions.', required{idx});
                end
            end

            vectorSize = [];
            for idx = 1:numel(required)
                value = params.(required{idx});
                validateattributes(value, {'numeric'}, {'real', 'nonnegative'});

                if ~isscalar(value)
                    thisSize = size(value);
                    if isempty(vectorSize)
                        vectorSize = thisSize;
                    elseif ~isequal(thisSize, vectorSize)
                        error('AnalyticalCylinder3D:ShapeParameters:LengthMismatch', ...
                            'Vector-valued radius/length must share the same size.');
                    end
                end
            end
        end

        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            params = obj.getShapeParameters();
            kr = sqrt(kx.^2 + ky.^2);
            arg = 2*pi*params.radius_mm .* kr;

            radial = (params.radius_mm .* besselj(1, arg)) ./ kr;
            if(isscalar(params.radius_mm))
                radial(kr == 0) = pi .* params.radius_mm.^2;
            else
                radial(kr == 0) = pi .* params.radius_mm(kr == 0).^2;
            end
            axial = params.length_mm .* sinc(params.length_mm .* kz);

            S_body = radial .* axial;
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            params = obj.getShapeParameters();
            inRadius = (xb.^2 + yb.^2) <= (params.radius_mm.^2);
            inHeight = abs(zb) <= (params.length_mm ./ 2);
            percent = double(inRadius & inHeight);
        end

        function params = defaultCylinderParameters(~)
            params = struct('radius_mm', 1, 'length_mm', 1);
        end
    end
end
