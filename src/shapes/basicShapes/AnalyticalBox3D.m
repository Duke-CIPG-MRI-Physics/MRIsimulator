classdef AnalyticalBox3D < AnalyticalShape3D
    % AnalyticalBox3D
    %   Rectangular prism aligned with BODY axes. Dimensions can be provided as
    %   static scalars/vectors, a struct of those values, or as a function
    %   handle returning such a struct for time-varying shapes.

    methods
        function obj = AnalyticalBox3D(intensity, shapeParameters)
            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = obj.defaultBoxParameters();
            end
            if nargin < 1 || isempty(intensity)
                intensity = 1;
            end

            shapeParameters = AnalyticalShape3D.ensurePoseFields(shapeParameters);
            obj@AnalyticalShape3D(intensity, shapeParameters);
        end
    end

    methods (Access = protected)
        function validateParameterFields(obj, params)
            validateParameterFields@AnalyticalShape3D(obj, params);

            required = {'Lx_mm', 'Ly_mm', 'Lz_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalBox3D:ShapeParameters:MissingField', ...
                        'Field %s is required for box dimensions.', required{idx});
                end

                value = params.(required{idx});
                if isnumeric(value)
                    validateattributes(value, {'numeric'}, {'real', 'nonnegative', 'finite'});
                elseif ~isa(value, 'function_handle')
                    error('AnalyticalBox3D:ShapeParameters:InvalidType', ...
                        'Dimensions must be numeric or function handles.');
                end
            end
        end

        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            params = obj.getShapeParameters();
            S_body = params.Lx_mm .* params.Ly_mm .* params.Lz_mm .* ... % volume
                sinc(kx .* params.Lx_mm) .* ...
                sinc(ky .* params.Ly_mm) .* ...
                sinc(kz .* params.Lz_mm);
        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            params = obj.getShapeParameters();
            inX = abs(xb) <= params.Lx_mm ./ 2;
            inY = abs(yb) <= params.Ly_mm ./ 2;
            inZ = abs(zb) <= params.Lz_mm ./ 2;
            percent = double(inX & inY & inZ);
        end

        function params = defaultBoxParameters(~)
            params = struct('Lx_mm', 1, 'Ly_mm', 1, 'Lz_mm', 1);
        end
    end
end
