classdef MultipleMaterialPhantom < AnalyticalShape3D
    % MultipleMaterialPhantom
    %   Container for a collection of AnalyticalShape3D objects representing
    %   different materials. Evaluates k-space by summing the contributions
    %   from each contained shape (each shape carries its own intensity).

    properties (Access = protected)
        shapes (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
    end

    methods
        function obj = MultipleMaterialPhantom(shapes)
            obj@AnalyticalShape3D(1, [0 0 0], [0 0 0]);

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

    methods (Access = protected, Hidden = false)
        function S = bodyKspace(obj, kx_body, ky_body, kz_body)
            % bodyKspace  BODY-frame analytic FT (no intensity scaling).
            %   Implemented as the sum of each contained shape's WORLD-frame
            %   k-space because the phantom itself does not apply any
            %   additional transform beyond the individual shapes.

            if isempty(obj.shapes)
                S = zeros(size(kx_body));
                return;
            end

            S = zeros(size(kx_body));
            for idx = 1:numel(obj.shapes)
                S = S + obj.shapes(idx).kspace(kx_body, ky_body, kz_body);
            end
        end

        function percent = percentInsideShape(obj, xb, yb, zb)
            % percentInsideShape  Fraction of each voxel occupied by the
            % composite phantom (0–1). Overlaps are handled by computing the
            % complement of the product of "outside" probabilities for each
            % contained shape.

            if isempty(obj.shapes)
                percent = zeros(size(xb));
                return;
            end

            percent = zeros(size(xb));
            for idx = 1:numel(obj.shapes)
                pctShape = obj.shapes(idx).percentInsideShape(xb, yb, zb);
                percent = 1 - (1 - percent) .* (1 - pctShape);
            end
        end
    end
end
