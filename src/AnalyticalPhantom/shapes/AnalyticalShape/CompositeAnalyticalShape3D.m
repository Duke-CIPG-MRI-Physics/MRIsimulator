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
    %       CompositeAnalyticalShape3D instance containing that handle. The
    %       composite listens for component shapeChanged events and forwards a
    %       shapeChanged notification when any component changes.
    %
    %   Constructor:
    %       obj = CompositeAnalyticalShape3D(additiveComponents, subtractiveComponents, ...
    %               intensity, center, rollPitchYaw)
    %
    %   Methods:
    %       addComponent(shape)      % append to additiveComponents
    %       subtractComponent(shape) % append to subtractiveComponents

    properties (Access = protected)
        additiveComponents (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
        subtractiveComponents (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
    end

    properties (Access = private)
        componentListeners = event.listener.empty(1,0);
    end

    methods
        function obj = CompositeAnalyticalShape3D(additiveComponents, subtractiveComponents, ...
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
            obj.registerComponentListeners(shape);
            obj.markShapeChanged();
        end

        function subtractComponent(obj, shape)
            arguments
                obj
                shape (1,:) AnalyticalShape3D
            end

            obj.subtractiveComponents = [obj.subtractiveComponents, shape];
            obj.registerComponentListeners(shape);
            obj.markShapeChanged();
        end
    end

    methods (Access = private)
        function registerComponentListeners(obj, components)
            for idx = 1:numel(components)
                comp = components(idx);
                lh = comp.addShapeChangedListener(@(~,~) obj.markShapeChanged()); %#ok<AGROW>
                obj.componentListeners(end+1) = lh;
            end
        end
    end

    methods (Access = protected)
        function S = bodyKspace(obj, kx_body, ky_body, kz_body)
            arguments
                obj
                kx_body double
                ky_body double
                kz_body double
            end

            if ~isequal(size(kx_body), size(ky_body), size(kz_body))
                error('CompositeAnalyticalShape3D:bodyKspace:SizeMismatch', ...
                    'kx_body, ky_body, kz_body must have identical sizes.');
            end

            
            S_body = zeros(size(kx_body));
            for idx = 1:numel(obj.additiveComponents)
                S_body = S_body + obj.additiveComponents(idx).kspace_shapeOnly(kx_body, ky_body, kz_body);
            end

            for idx = 1:numel(obj.subtractiveComponents)
                S_body = S_body - obj.subtractiveComponents(idx).kspace_shapeOnly(kx_body, ky_body, kz_body);
            end

            % WORLD translation phase
            c = obj.getCenter();
            if any(c ~= 0, 'all')
                if size(c, 2) ~= 3
                    error('AnalyticalShape3D:kspace_shapeOnly:CenterSizeMismatch', ...
                        'Center must have 3 columns for x, y, z.');
                end

                cx = obj.requireScalarOrSize(c(:,1), kx_body, 'centerX');
                cy = obj.requireScalarOrSize(c(:,2), ky_body, 'centerY');
                cz = obj.requireScalarOrSize(c(:,3), kz_body, 'centerZ');

                phase = exp(-1i * 2*pi * ( ...
                        kx_body .* cx + ...
                        ky_body .* cy + ...
                        kz_body .* cz));
                S = S_body .* phase;
            else
                S = S_body;
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
    end
end
