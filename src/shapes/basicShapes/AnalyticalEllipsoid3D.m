classdef AnalyticalEllipsoid3D < AnalyticalShape3D
    % AnalyticalEllipsoid3D
    %   Solid ellipsoid aligned with BODY axes. Semi-axes can be static
    %   scalars/vectors, a struct of those values, or supplied via a function
    %   handle returning such a struct for time-varying shapes.

    methods
        function obj = AnalyticalEllipsoid3D(intensity, shapeParameters)
            if nargin < 1 || isempty(intensity)
                intensity = 1;
            end

            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = AnalyticalEllipsoid3D.defaultEllipsoidParameters();
            end

            obj@AnalyticalShape3D(intensity, shapeParameters);
        end
    end

    methods (Access = protected)
        function validateParameterFields(obj, params)
            validateParameterFields@AnalyticalShape3D(obj, params);

            required = {'a_mm', 'b_mm', 'c_mm'};
            for idx = 1:numel(required)
                if ~isfield(params, required{idx})
                    error('AnalyticalEllipsoid3D:ShapeParameters:MissingField', ...
                        'Field %s is required for ellipsoid dimensions.', required{idx});
                end

                value = params.(required{idx});
                if isnumeric(value)
                    validateattributes(value, {'numeric'}, {'real', 'nonnegative', 'finite'});
                elseif ~isa(value, 'function_handle')
                    error('AnalyticalEllipsoid3D:ShapeParameters:InvalidType', ...
                        'Dimensions must be numeric or function handles.');
                end
            end
        end

        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            params = obj.getAxesParameters();
            kScaled = sqrt((params.a_mm .* kx).^2 + (params.b_mm .* ky).^2 + (params.c_mm .* kz).^2);
            arg = 2 * pi .* kScaled;

            % Spherical Bessel j1(x) = (sin(x) - x cos(x)) / x^2
            numerator = sin(arg) - arg .* cos(arg);
            denom = (arg).^3;

            vol = (4/3) * pi .* params.a_mm .* params.b_mm .* params.c_mm;
            S_body = vol .* 3 .* numerator ./ denom;
            if(isscalar(vol))
                S_body(kScaled == 0) = vol;
            else
                S_body(kScaled == 0) = vol(kScaled == 0);
            end

        end

        function percent = percentInsideBody(obj, xb, yb, zb)
            params = obj.getAxesParameters();
            scaled = (xb ./ params.a_mm).^2 + (yb ./ params.b_mm).^2 + (zb ./ params.c_mm).^2;
            percent = double(scaled <= 1);
        end

        function params = getAxesParameters(obj, varargin)
            params = obj.getShapeParameters(varargin{:});
        end
    end

    methods (Static, Access = protected)
        function params = defaultEllipsoidParameters()
            params = struct('a_mm', 1, 'b_mm', 1, 'c_mm', 1);
        end
    end
end
