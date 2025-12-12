classdef MultipleMaterialPhantom < AnalyticalShape3D
    % MultipleMaterialPhantom
    %   Container for a collection of AnalyticalShape3D objects representing
    %   different materials. Evaluates k-space by summing the contributions
    %   from each contained shape (each shape carries its own intensity).
    %
    %   Constructor:
    %       phantom = MultipleMaterialPhantom(shapes, shapeParameters)

    properties (Access = protected)
        shapes  % No type constraint - validated in setShapes() method
    end

    methods
        function obj = MultipleMaterialPhantom(shapes, shapeParameters)
            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = AnalyticalShape3D.ensurePoseFields(struct());
            end
            if nargin < 1
                shapes = [];
            end

            obj@AnalyticalShape3D(1, shapeParameters);

            % Initialize shapes as empty array of the abstract class
            obj.shapes = AnalyticalShape3D.empty(1,0);

            if ~isempty(shapes)
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

            params = obj.getShapeParameters();
            pose = obj.extractPose(params, kx);
            rpy = [pose.roll_deg, pose.pitch_deg, pose.yaw_deg];
            noRotation = ~any(rpy(:));

            if noRotation
                % Fast path: BODY and WORLD frames align
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                % WORLD → BODY transform
                [kxb, kyb, kzb] = obj.rotateWorldToBody(kx, ky, kz, pose.roll_deg, pose.pitch_deg, pose.yaw_deg);
            end

            S_body = zeros(size(kxb));
            for idx = 1:numel(obj.shapes)
                S_body = S_body + obj.shapes(idx).kspace(kxb, kyb, kzb);
            end

             % WORLD translation phase
            if any(pose.centerX(:) ~= 0) || any(pose.centerY(:) ~= 0) || any(pose.centerZ(:) ~= 0)
                phase = exp(-1i * 2*pi * ( ...
                        kx .* pose.centerX + ...
                        ky .* pose.centerY + ...
                        kz .* pose.centerZ));
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

            params = obj.getShapeParameters();
            pose = obj.extractPose(params, kx);
            rpy = [pose.roll_deg, pose.pitch_deg, pose.yaw_deg];
            noRotation = ~any(rpy(:));

            if noRotation
                % Fast path: BODY and WORLD frames align
                kxb = kx;
                kyb = ky;
                kzb = kz;
            else
                % WORLD → BODY transform
                [kxb, kyb, kzb] = obj.rotateWorldToBody(kx, ky, kz, pose.roll_deg, pose.pitch_deg, pose.yaw_deg);
            end

            S_body = zeros(size(kxb));
            for idx = 1:numel(obj.shapes)
                S_body = S_body + obj.shapes(idx).kspaceWorldPlacedShape(kxb, kyb, kzb);
            end

            % WORLD translation phase
            if any(pose.centerX(:) ~= 0) || any(pose.centerY(:) ~= 0) || any(pose.centerZ(:) ~= 0)
                phase = exp(-1i * 2*pi * ( ...
                        kx .* pose.centerX + ...
                        ky .* pose.centerY + ...
                        kz .* pose.centerZ));
                S = S_body .* phase;
            else
                S = S_body;
            end
        end
    end
    methods (Access = protected)
        function validateParameterFields(obj, params)
            validateParameterFields@AnalyticalShape3D(obj, params);
        end

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
                percent = zeros(size(xb));
                return;
            end

            percent = zeros(size(xb));
            for idx = 1:numel(obj.shapes)
                percent = percent + obj.shapes(idx).percentInsideBody(xb, yb, zb);
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
