classdef MultiIntensityShapeGroup3D < AnalyticalShape3D
    % MultiIntensityShapeGroup3D
    %   Groups AnalyticalShape3D objects (or SharedIntensityShapeGroup3D
    %   instances) that each retain their own intensity. K-space is the sum of
    %   all children, optionally sharing a common WORLD center/orientation.
    %
    %   Differences vs SharedIntensityShapeGroup3D:
    %     • MultiIntensityShapeGroup3D preserves each child's intensity and is
    %       suited for phantoms with multiple labeled materials.
    %     • SharedIntensityShapeGroup3D ignores child intensities and applies a
    %       single scaling to the boolean blend of additive/subtractive parts.
    %     • Both containers allow grouped pose so nested shapes move together.

    properties (Access = protected)
        shapes (1,:) AnalyticalShape3D = AnalyticalShape3D.empty;
    end

    methods
        function obj = MultiIntensityShapeGroup3D(shapes, center, rollPitchYaw)
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

        function S = kspaceWorldGeometryScaled(obj, kx, ky, kz)
            % kspaceWorldGeometryScaled  Sum analytic k-space over all materials,
            %   including each shape's own intensity scaling.
            %   S = kspaceWorldGeometryScaled(obj, kx, ky, kz)
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
                S_body = S_body + obj.shapes(idx).kspaceWorldGeometryScaled(kxb, kyb, kzb);
            end

             % WORLD translation phase
            c = obj.getCenter();
            if any(c(:) ~= 0)
                if size(c, 2) ~= 3
                    error('MultiIntensityShapeGroup3D:kspaceWorldGeometryScaled:CenterSizeMismatch', ...
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

        function S = kspaceWorldGeometry(obj, kx, ky, kz)
            % kspaceWorldGeometry  Sum geometry-only k-space for all materials.
            %   S = kspaceWorldGeometry(obj, kx, ky, kz)
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
                S_body = S_body + obj.shapes(idx).kspaceWorldGeometry(kxb, kyb, kzb);
            end

            % WORLD translation phase
            c = obj.getCenter();
            if any(c(:) ~= 0)
                if size(c, 2) ~= 3
                    error('MultiIntensityShapeGroup3D:kspaceWorldGeometry:CenterSizeMismatch', ...
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
        function S = kspaceBodyGeometry(obj, kx_body, ky_body, kz_body)
            % kspaceBodyGeometry  BODY-frame analytic FT (no intensity scaling).
            %   Implemented as the sum of each contained shape's WORLD-frame
            %   k-space because the phantom itself does not apply any
            %   additional transform beyond the individual shapes.

            if isempty(obj.shapes)
                S = zeros(size(kx_body));
                return;
            end

            S = zeros(size(kx_body));
            for idx = 1:numel(obj.shapes)
                S = S + obj.shapes(idx).kspaceWorldGeometryScaled(kx_body, ky_body, kz_body);
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
