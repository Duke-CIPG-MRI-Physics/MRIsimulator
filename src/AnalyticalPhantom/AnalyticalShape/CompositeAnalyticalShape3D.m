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
    %       kspace_shapeOnly(kx, ky, kz) % sum additive minus subtractive FTs

    properties (Access = protected)
        additiveComponents (1,:) AnalyticalShape3D = AnalyticalShape3D.empty(1,0);
        subtractiveComponents (1,:) AnalyticalShape3D = AnalyticalShape3D.empty(1,0);
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

    methods
        function S = kspace_shapeOnly(obj, kx, ky, kz)
            arguments
                obj
                kx double
                ky double
                kz double
            end

            if ~isequal(size(kx), size(ky), size(kz))
                error('CompositeAnalyticalShape3D:kspace_shapeOnly:SizeMismatch', ...
                    'kx, ky, kz must have identical sizes.');
            end

            S = zeros(size(kx));

            for idx = 1:numel(obj.additiveComponents)
                S = S + obj.additiveComponents(idx).kspace_shapeOnly(kx, ky, kz);
            end

            for idx = 1:numel(obj.subtractiveComponents)
                S = S - obj.subtractiveComponents(idx).kspace_shapeOnly(kx, ky, kz);
            end
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
end
