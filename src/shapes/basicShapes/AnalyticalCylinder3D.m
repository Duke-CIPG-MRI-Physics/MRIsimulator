classdef AnalyticalCylinder3D < AnalyticalShape3D
    % AnalyticalCylinder3D
    %   Finite cylinder aligned with the BODY +z axis. Radius/length can be
    %   static scalars/vectors, a struct of those values, or supplied by a
    %   function handle returning such a struct for time-varying geometries.

    methods
        function obj = AnalyticalCylinder3D(intensity, shapeParameters)
            if nargin < 1 || isempty(intensity)
                intensity = 1;
            end

            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = AnalyticalCylinder3D.defaultCylinderParameters();
            end

            obj@AnalyticalShape3D(intensity, shapeParameters);
        end
    end

    methods (Access = protected)
        function validateParameterFields(obj, params)
            validateParameterFields@AnalyticalShape3D(obj, params);

            required = {'radius_mm', 'length_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalCylinder3D:ShapeParameters:MissingField', ...
                        'Field %s is required for cylinder dimensions.', required{idx});
                end

                value = params.(required{idx});
                if isnumeric(value)
                    validateattributes(value, {'numeric'}, {'real', 'nonnegative', 'finite'});
                elseif ~isa(value, 'function_handle')
                    error('AnalyticalCylinder3D:ShapeParameters:InvalidType', ...
                        'Dimensions must be numeric or function handles.');
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
    end

    methods (Static, Access = protected)
        function params = defaultCylinderParameters()
            params = struct('radius_mm', 1, 'length_mm', 1);
        end
    end
end
