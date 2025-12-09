classdef MultipleMaterialPhantom < AnalyticalShape3D
    % MultipleMaterialPhantom
    %   Container for a collection of AnalyticalShape3D objects representing
    %   different materials. Evaluates k-space by summing the contributions
    %   from each contained shape (each shape carries its own intensity).

    properties (Access = protected)
        shapes (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
    end

    methods
        function obj = MultipleMaterialPhantom(shapes, center, rollPitchYaw)
            if nargin < 3 
                rollPitchYaw = [0, 0, 0];
            end
            if nargin < 2
                center = [0, 0, 0];
            end

            obj@AnalyticalShape3D(1, center, rollPitchYaw);

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
            % kspace  Sum placed-and-scaled k-space over all materials.
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

            rpy = obj.getRollPitchYaw();
            noRotation = ~any(rpy);

            if noRotation
                % Fast path: BODY and WORLD frames align
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                % WORLD → BODY transform
                inputSize = size(kx);
                K_world = [kx(:) ky(:) kz(:)].';
                R = obj.calculateRotationMatrix();
                K_body = R.' * K_world;

                kxb = reshape(K_body(1,:), inputSize);
                kyb = reshape(K_body(2,:), inputSize);
                kzb = reshape(K_body(3,:), inputSize);
            end

            S_body = zeros(size(kxb));
            for idx = 1:numel(obj.shapes)
                S_body = S_body + obj.shapes(idx).kspace(kxb, kyb, kzb);
            end

             % WORLD translation phase
            c = obj.getCenter();
            if any(c(:) ~= 0)
                if size(c, 2) ~= 3
                    error('AnalyticalShape3D:kspaceWorldPlacedShape:CenterSizeMismatch', ...
                        'Center must have 3 columns for x, y, z.');
                end

                cx = obj.requireScalarOrSize(c(:,1), kx, 'centerX');
                cy = obj.requireScalarOrSize(c(:,2), ky, 'centerY');
                cz = obj.requireScalarOrSize(c(:,3), kz, 'centerZ');

                phase = exp(-1i * 2*pi * ( ...
                        kx .* cx + ...
                        ky .* cy + ...
                        kz .* cz));
                S = S_body .* phase;
            else
                S = S_body;
            end
        end

        function S = kspaceWorldPlacedShape(obj, kx, ky, kz)
            % kspaceWorldPlacedShape  Sum geometry-only k-space for all materials.
            %   S = kspaceWorldPlacedShape(obj, kx, ky, kz)
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

            rpy = obj.getRollPitchYaw();
            noRotation = ~any(rpy);

            if noRotation
                % Fast path: BODY and WORLD frames align
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                % WORLD → BODY transform
                inputSize = size(kx);
                K_world = [kx(:) ky(:) kz(:)].';
                R = obj.calculateRotationMatrix();
                K_body = R.' * K_world;

                kxb = reshape(K_body(1,:), inputSize);
                kyb = reshape(K_body(2,:), inputSize);
                kzb = reshape(K_body(3,:), inputSize);
            end

            S_body = zeros(size(kxb));
            for idx = 1:numel(obj.shapes)
                S_body = S_body + obj.shapes(idx).kspaceWorldPlacedShape(kxb, kyb, kzb);
            end

            % WORLD translation phase
            c = obj.getCenter();
            if any(c(:) ~= 0)
                if size(c, 2) ~= 3
                    error('AnalyticalShape3D:kspaceWorldPlacedShape:CenterSizeMismatch', ...
                        'Center must have 3 columns for x, y, z.');
                end

                cx = obj.requireScalarOrSize(c(:,1), kx, 'centerX');
                cy = obj.requireScalarOrSize(c(:,2), ky, 'centerY');
                cz = obj.requireScalarOrSize(c(:,3), kz, 'centerZ');

                phase = exp(-1i * 2*pi * ( ...
                        kx .* cx + ...
                        ky .* cy + ...
                        kz .* cz));
                S = S_body .* phase;
            else
                S = S_body;
            end
        end
    end
    methods (Access = protected)
        function S = kspaceBaseShape(obj, kx_body, ky_body, kz_body)
            % kspaceBaseShape  BODY-frame analytic FT (no intensity scaling).
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

        
        function percent = percentInsideBody(obj, xb, yb, zb)
            if isempty(obj.shapes)
                S = zeros(size(kx_body));
                return;
            end

            S = zeros(size(kx_body));
            for idx = 1:numel(obj.shapes)
                S = S + obj.shapes(idx).percentInsideBody(kx, ky, kz);
            end
        end
  
    end

    methods
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
