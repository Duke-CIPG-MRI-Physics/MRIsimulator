classdef MultipleMaterialPhantom < handle
    % MultipleMaterialPhantom
    %   Container for a collection of AnalyticalShape3D objects representing
    %   different materials. Evaluates k-space by summing the contributions
    %   from each contained shape (each shape carries its own intensity).

    properties (Access = protected)
        shapes (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
    end

    methods
        function obj = MultipleMaterialPhantom(shapes)
            if nargin >= 1 && ~isempty(shapes)
                obj.setShapes(shapes);
            end
        end

        function setShapes(obj, shapes)
            arguments
                obj
                shapes (1,:) AnalyticalShape3D
            end

            obj.shapes = shapes;
        end

        function addShape(obj, shape)
            arguments
                obj
                shape (1,:) AnalyticalShape3D
            end

            obj.shapes = [obj.shapes, shape];
        end

        function S = kspace(obj, kx, ky, kz)
            % kspace  Sum analytic k-space over all materials.
            %   S = kspace(obj, kx, ky, kz)
            %
            %   Inputs:
            %       kx, ky, kz : WORLD frequencies [cycles/mm], same size.
            %
            %   Output:
            %       S : complex double, same size as inputs.

            if isempty(obj.shapes)
                S = zeros(size(kx));
                return;
            end

            S = zeros(size(kx));
            for idx = 1:numel(obj.shapes)
                S = S + obj.shapes(idx).kspace(kx, ky, kz);
            end
        end

        function S = kspace_shapeOnly(obj, kx, ky, kz)
            % kspace_shapeOnly  Sum geometry-only k-space for all materials.
            %   S = kspace_shapeOnly(obj, kx, ky, kz)
            %
            %   Inputs:
            %       kx, ky, kz : WORLD frequencies [cycles/mm], same size.
            %
            %   Output:
            %       S : complex double, same size as inputs.

            if isempty(obj.shapes)
                S = zeros(size(kx));
                return;
            end

            S = zeros(size(kx));
            for idx = 1:numel(obj.shapes)
                S = S + obj.shapes(idx).kspace_shapeOnly(kx, ky, kz);
            end
        end
    end
end
