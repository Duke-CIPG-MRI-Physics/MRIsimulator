classdef AnalyticalEllipticalCylinder3D < AnalyticalShape3D
    % AnalyticalEllipticalCylinder3D
    %   Finite cylinder with an elliptical cross section aligned to BODY
    %   axes. Semi-axes/length can be static scalars/vectors, provided as a
    %   struct, or computed via a function handle returning such a struct for
    %   time-varying geometry.

    methods
        function obj = AnalyticalEllipticalCylinder3D(intensity, shapeParameters)
            if nargin < 1 || isempty(intensity)
                intensity = 1;
            end

            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = AnalyticalEllipticalCylinder3D.defaultEllipticalCylinderParameters();
            end

            obj@AnalyticalShape3D(intensity, shapeParameters);
        end
    end

    methods (Access = protected)
        function validateParameterFields(obj, params)
            validateParameterFields@AnalyticalShape3D(obj, params);

            required = {'a_mm', 'b_mm', 'length_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalEllipticalCylinder3D:ShapeParameters:MissingField', ...
                        'Field %s is required for elliptical cylinder dimensions.', required{idx});
                end

                value = params.(required{idx});
                if isnumeric(value)
                    validateattributes(value, {'numeric'}, {'real', 'nonnegative', 'finite'});
                elseif ~isa(value, 'function_handle')
                    error('AnalyticalEllipticalCylinder3D:ShapeParameters:InvalidType', ...
                        'Dimensions must be numeric or function handles.');
                end
            end
        end

        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            params = obj.getShapeParameters();
            kr = sqrt((params.a_mm .* kx).^2 + (params.b_mm .* ky).^2);
            arg = 2 * pi .* kr;

            radial = (params.a_mm .* params.b_mm .* besselj(1, arg)) ./ kr;
            if(isscalar(params.a_mm) & isscalar(params.b_mm))
                radial(kr == 0) = pi .* params.a_mm .* params.b_mm;
            else
                radial(kr == 0) = pi .* params.a_mm(kr == 0) .* params.b_mm(kr == 0);
            end

            axial = params.length_mm .* sinc(params.length_mm .* kz);

            S_body = radial .* axial;
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            params = obj.getShapeParameters();
            inEllipse = (xb ./ params.a_mm).^2 + (yb ./ params.b_mm).^2 <= 1;
            inHeight = abs(zb) <= params.length_mm ./ 2;
            percent = double(inEllipse & inHeight);
        end
    end

    methods (Static, Access = protected)
        function params = defaultEllipticalCylinderParameters()
            params = struct('a_mm', 1, 'b_mm', 1, 'length_mm', 1);
        end
    end
end
