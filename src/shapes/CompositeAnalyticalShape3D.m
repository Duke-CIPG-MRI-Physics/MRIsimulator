classdef CompositeAnalyticalShape3D < AnalyticalShape3D
    % CompositeAnalyticalShape3D
    %   Aggregates additive and subtractive AnalyticalShape3D components.
    %
    %   Notes:
    %     • Component intensities are ignored; this class treats components as
    %       geometry-only building blocks and applies its own shapeIntensity to
    %       the combined result.
    %     • Components are handle objects. Any external mutation (geometry or
    %       intensity) to a shared component will affect every
    %       CompositeAnalyticalShape3D instance containing that handle.
    %
    %   Constructor:
    %       obj = CompositeAnalyticalShape3D(additiveComponents, subtractiveComponents, ...
    %               intensity, shapeParameters)
    %
    %   Methods:
    %       addComponent(shape)      % append to additiveComponents
    %       subtractComponent(shape) % append to subtractiveComponents

    properties (Access = protected)
        additiveComponents  % No type constraint
        subtractiveComponents  % No type constraint
    end

    methods
        function obj = CompositeAnalyticalShape3D(additiveComponents, subtractiveComponents, ...
                intensity, shapeParameters)
            if nargin < 1 || isempty(additiveComponents)
                additiveComponents = AnalyticalShape3D.empty(1,0);
            end

            if nargin < 2 || isempty(subtractiveComponents)
                subtractiveComponents = AnalyticalShape3D.empty(1,0);
            end

            if nargin < 3 || isempty(intensity)
                intensity = 1;
            end

            if nargin < 4 || isempty(shapeParameters)
                shapeParameters = AnalyticalShape3D.ensurePoseFields(struct());
            else
                shapeParameters = AnalyticalShape3D.ensurePoseFields(shapeParameters);
            end

            obj@AnalyticalShape3D(intensity, shapeParameters);

            % Initialize component arrays as empty
            obj.additiveComponents = AnalyticalShape3D.empty(1,0);
            obj.subtractiveComponents = AnalyticalShape3D.empty(1,0);

            if nargin >= 1 && ~isempty(additiveComponents)
                obj.addComponent(additiveComponents);
            end

            if nargin >= 2 && ~isempty(subtractiveComponents)
                obj.subtractComponent(subtractiveComponents);
            end
        end

        function addComponent(obj, shape)
            arguments
                obj
                shape (1,:) AnalyticalShape3D
            end

            obj.additiveComponents = [obj.additiveComponents, shape];
        end

        function subtractComponent(obj, shape)
            arguments
                obj
                shape (1,:) AnalyticalShape3D
            end

            obj.subtractiveComponents = [obj.subtractiveComponents, shape];
        end
    end

    methods
        function image = estimateImage(obj, xMesh, yMesh, zMesh)
            arguments
                obj
                xMesh double
                yMesh double
                zMesh double
            end

            if ~isequal(size(xMesh), size(yMesh), size(zMesh))
                error('CompositeAnalyticalShape3D:estimateImage:SizeMismatch', ...
                    'xMesh, yMesh, zMesh must have identical sizes.');
            end

            image = zeros(size(xMesh));

            for idx = 1:numel(obj.additiveComponents)
                image = image + obj.additiveComponents(idx).estimateImage(xMesh, yMesh, zMesh);
            end

            for idx = 1:numel(obj.subtractiveComponents)
                image = image - obj.subtractiveComponents(idx).estimateImage(xMesh, yMesh, zMesh);
            end
        end

        function mask = percentInsideShape(obj, xb, yb, zb)
            additiveMask = false(size(xb));
            for idx = 1:numel(obj.additiveComponents)
                additiveMask = additiveMask | (obj.additiveComponents(idx).percentInsideShape(xb, yb, zb) ~= 0);
            end

            subtractiveMask = false(size(xb));
            for idx = 1:numel(obj.subtractiveComponents)
                subtractiveMask = subtractiveMask | (obj.subtractiveComponents(idx).percentInsideShape(xb, yb, zb) ~= 0);
            end

            mask = additiveMask & ~subtractiveMask;
        end
    end

    methods (Access = protected)
        function validateParameterFields(obj, params)
            validateParameterFields@AnalyticalShape3D(obj, params);
        end

        function S_body = kspaceBaseShape(obj, kx, ky, kz)
            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('CompositeAnalyticalShape3D:kspaceBaseShape:SizeMismatch', ...
                    'kx, ky, kz must have identical sizes.');
            end

            S_body = zeros(size(kx));
            for idx = 1:numel(obj.additiveComponents)
                S_body = S_body + obj.additiveComponents(idx).kspaceWorldPlacedShape(kx, ky, kz);
            end

            for idx = 1:numel(obj.subtractiveComponents)
                S_body = S_body - obj.subtractiveComponents(idx).kspaceWorldPlacedShape(kx, ky, kz);
            end
        end

    end
end
