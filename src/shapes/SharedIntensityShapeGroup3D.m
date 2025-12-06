classdef SharedIntensityShapeGroup3D < AnalyticalShape3D
    % SharedIntensityShapeGroup3D
    %   Groups shapes into additive/subtractive Boolean blends that share a
    %   common WORLD pose and a single intensity scale.
    %
    %   Differences vs MultiIntensityShapeGroup3D:
    %     • SharedIntensityShapeGroup3D ignores component intensities and
    %       treats all children as geometry-only building blocks (the group's
    %       own intensity scales the final result).
    %     • MultiIntensityShapeGroup3D keeps each child's own intensity and
    %       sums them directly; use that when you want multiple labeled
    %       materials.
    %     • Both containers allow grouped center/orientation so nested shapes
    %       move together.
    %
    %   Notes:
    %     • Components are handle objects. Any external mutation (geometry or
    %       intensity) to a shared component will affect every
    %       SharedIntensityShapeGroup3D instance containing that handle. The
    %       composite queries component geometry directly when evaluated rather
    %       than relying on listener callbacks.
    %
    %   Constructor:
    %       obj = SharedIntensityShapeGroup3D(additiveComponents, subtractiveComponents, ...
    %               intensity, center, rollPitchYaw)
    %
    %   Methods:
    %       addComponent(shape)      % append to additiveComponents
    %       subtractComponent(shape) % append to subtractiveComponents

    properties (Access = protected)
        additiveComponents (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
        subtractiveComponents (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
    end

    methods
        function obj = SharedIntensityShapeGroup3D(additiveComponents, subtractiveComponents, ...
                intensity, center, rollPitchYaw)
            obj@AnalyticalShape3D(intensity, center, rollPitchYaw);

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
            obj.markShapeChanged();
        end

        function subtractComponent(obj, shape)
            arguments
                obj
                shape (1,:) AnalyticalShape3D
            end

            obj.subtractiveComponents = [obj.subtractiveComponents, shape];
            obj.markShapeChanged();
        end
    end

    methods (Access = protected)
        function S_body = kspaceBodyGeometry(obj, kx, ky, kz)
            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('SharedIntensityShapeGroup3D:kspaceBodyGeometry:SizeMismatch', ...
                    'kx, ky, kz must have identical sizes.');
            end

            S_body = zeros(size(kx));
            for idx = 1:numel(obj.additiveComponents)
                S_body = S_body + obj.additiveComponents(idx).kspaceWorldGeometry(kx, ky, kz);
            end

            for idx = 1:numel(obj.subtractiveComponents)
                S_body = S_body - obj.subtractiveComponents(idx).kspaceWorldGeometry(kx, ky, kz);
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

    methods
        function image = estimateImage(obj, xMesh, yMesh, zMesh)
            arguments
                obj
                xMesh double
                yMesh double
                zMesh double
            end

            if ~isequal(size(xMesh), size(yMesh), size(zMesh))
                error('SharedIntensityShapeGroup3D:estimateImage:SizeMismatch', ...
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
    end
end
