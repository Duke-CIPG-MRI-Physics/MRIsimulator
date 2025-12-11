classdef AnalyticalBox3D < AnalyticalShape3D
    % AnalyticalBox3D
    %   Rectangular prism aligned with BODY axes. Dimensions can be provided as
    %   static scalars/vectors, a struct of those values, or as a function
    %   handle returning such a struct for time-varying shapes.

    methods
        function obj = AnalyticalBox3D(shapeParameters, intensity, center, rollPitchYaw)
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
                shapeParameters = obj.defaultBoxParameters();
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
                error('AnalyticalBox3D:ShapeParameters:InvalidType', ...
                    'Shape parameters must be provided as a struct.');
            end

            required = {'Lx_mm', 'Ly_mm', 'Lz_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalBox3D:ShapeParameters:MissingField', ...
                        'Field %s is required for box dimensions.', required{idx});
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
                        error('AnalyticalBox3D:ShapeParameters:LengthMismatch', ...
                            'All vector-valued dimensions must share the same size.');
                    end
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
