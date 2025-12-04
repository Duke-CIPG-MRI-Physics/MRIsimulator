classdef CompositeAnalyticalShape3D < AnalyticalShape3D
    % CompositeAnalyticalShape3D
    %   Aggregates additive and subtractive AnalyticalShape3D components.
    %
    %   Notes:
    %     • Component intensities are ignored; this class treats components as
    %       geometry-only building blocks and applies its own shapeIntensity to
    %       the combined result.
    %     • Component centers and orientations are interpreted relative to the
    %       composite's BODY frame so that nesting preserves pose hierarchies.
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

        function componentKspace = evaluateComponentKspace(~, component, kx_body, ky_body, kz_body)
            rpy = component.getRollPitchYaw();
            noRotation = ~any(rpy);

            if noRotation
                kx_componentBody = kx_body;
                ky_componentBody = ky_body;
                kz_componentBody = kz_body;
            else
                inputSize = size(kx_body);
                K_body = [kx_body(:) ky_body(:) kz_body(:)].';
                R_component = component.calculateRotationMatrix();
                K_componentBody = R_component.' * K_body;

                kx_componentBody = reshape(K_componentBody(1,:), inputSize);
                ky_componentBody = reshape(K_componentBody(2,:), inputSize);
                kz_componentBody = reshape(K_componentBody(3,:), inputSize);
            end

            componentBodyKspace = component.bodyKspace(kx_componentBody, ky_componentBody, kz_componentBody);

            componentCenter = component.getCenter();
            if size(componentCenter, 2) ~= 3
                error('CompositeAnalyticalShape3D:bodyKspace:CenterSizeMismatch', ...
                    'Component center must have 3 columns for x, y, z.');
            end

            cx = component.requireScalarOrSize(componentCenter(:,1), kx_body, 'componentCenterX');
            cy = component.requireScalarOrSize(componentCenter(:,2), ky_body, 'componentCenterY');
            cz = component.requireScalarOrSize(componentCenter(:,3), kz_body, 'componentCenterZ');

            phase = exp(-1i * 2*pi * ( ...
                    kx_body .* cx + ...
                    ky_body .* cy + ...
                    kz_body .* cz));

            componentKspace = componentBodyKspace .* phase;
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

            S = zeros(size(kx_body));

            for idx = 1:numel(obj.additiveComponents)
                component = obj.additiveComponents(idx);
                componentKspace = obj.evaluateComponentKspace(component, kx_body, ky_body, kz_body);
                S = S + componentKspace;
            end

            for idx = 1:numel(obj.subtractiveComponents)
                component = obj.subtractiveComponents(idx);
                componentKspace = obj.evaluateComponentKspace(component, kx_body, ky_body, kz_body);
                S = S - componentKspace;
            end
        end

        function mask = percentInsideShape(obj, xb, yb, zb)
            additiveMask = false(size(xb));
            for idx = 1:numel(obj.additiveComponents)
                component = obj.additiveComponents(idx);
                [xb_comp, yb_comp, zb_comp] = component.worldToBodyPoints(xb, yb, zb);
                additiveMask = additiveMask | (component.percentInsideShape(xb_comp, yb_comp, zb_comp) ~= 0);
            end

            subtractiveMask = false(size(xb));
            for idx = 1:numel(obj.subtractiveComponents)
                component = obj.subtractiveComponents(idx);
                [xb_comp, yb_comp, zb_comp] = component.worldToBodyPoints(xb, yb, zb);
                subtractiveMask = subtractiveMask | (component.percentInsideShape(xb_comp, yb_comp, zb_comp) ~= 0);
            end

            mask = additiveMask & ~subtractiveMask;
        end
    end

end
