classdef AnalyticalEllipticalCylinder3D < AnalyticalShape3D
    % AnalyticalEllipticalCylinder3D
    %   Finite cylinder with an elliptical cross section aligned to BODY
    %   axes. Semi-axes/length can be static scalars/vectors, provided as a
    %   struct, or computed via a function handle returning such a struct for
    %   time-varying geometry.

    methods
        function obj = AnalyticalEllipticalCylinder3D(shapeParameters, intensity, center, rollPitchYaw)
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
                shapeParameters = obj.defaultEllipticalCylinderParameters();
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
                error('AnalyticalEllipticalCylinder3D:ShapeParameters:InvalidType', ...
                    'Shape parameters must be provided as a struct.');
            end

            required = {'a_mm', 'b_mm', 'length_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalEllipticalCylinder3D:ShapeParameters:MissingField', ...
                        'Field %s is required for elliptical cylinder dimensions.', required{idx});
                end
            end

            vectorSize = [];
            for idx = 1:numel(required)
                value = params.(required{idx});
                validateattributes(value, {'numeric'}, {'real', 'nonnegative'});
                if ~(isscalar(value) || isvector(value))
                    error('AnalyticalEllipticalCylinder3D:ShapeParameters:InvalidShape', ...
                        '%s must be scalar or vector-valued.', required{idx});
                end

                if ~isscalar(value)
                    thisSize = size(value);
                    if isempty(vectorSize)
                        vectorSize = thisSize;
                    elseif ~isequal(thisSize, vectorSize)
                        error('AnalyticalEllipticalCylinder3D:ShapeParameters:LengthMismatch', ...
                            'Vector-valued semi-axes and length must share the same size.');
                    end
                end
            end
        end

        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            params = obj.getShapeParameters();
            kr = sqrt((params.a_mm .* kx).^2 + (params.b_mm .* ky).^2);
            arg = 2 * pi .* kr;

            radial = (params.a_mm .* params.b_mm .* besselj(1, arg)) ./ kr;
            radial(kr == 0) = pi .* params.a_mm .* params.b_mm;

            axial = params.length_mm .* sinc(params.length_mm .* kz);

            S_body = radial .* axial;
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            params = obj.getShapeParameters();
            inEllipse = (xb ./ params.a_mm).^2 + (yb ./ params.b_mm).^2 <= 1;
            inHeight = abs(zb) <= params.length_mm ./ 2;
            percent = double(inEllipse & inHeight);
        end

        function params = defaultEllipticalCylinderParameters(~)
            params = struct('a_mm', 1, 'b_mm', 1, 'length_mm', 1);
        end
    end
end
