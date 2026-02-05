classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   single lesion in the right breast. Geometry and intensities follow
    %   the original demo_analyticalBreastPhantom.m setup. Heart geometry leverages
    %   cardiac_ellipsoid_waveform, which evaluates long time vectors in
    %   ~0.001 Gb, 15,625-sample chunks to keep memory usage minimal.

    properties (Access = protected)
        phantomParams
        beatingHeart
        rightLung
        leftLung
        fatInner
        fatOuter
        fatComposite
        tissueComposite
        thorax
        breastRight
        breastLeft
        leftAndRightBreastTissue
        lesionRight
    end

    methods
        function obj = BreastPhantom(params)
            arguments
                params struct = createBreastPhantomParams();
            end

            obj@MultipleMaterialPhantom();
            params = BreastPhantom.normalizeLesionParams(params);
            obj.phantomParams = params;

            % Convenience arrays
            centered = [0, 0, 0];
            notRotated = [0, 0, 0];

            placeholderEllipsoid = BreastPhantom.addPoseToParameters( ...
                struct('a_mm', 1, 'b_mm', 1, 'c_mm', 1), ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.beatingHeart = AnalyticalEllipsoid3D(params.heartIntensity, placeholderEllipsoid);
            obj.rightLung = AnalyticalEllipsoid3D(params.lungIntensity, placeholderEllipsoid);
            obj.leftLung = AnalyticalEllipsoid3D(params.lungIntensity, placeholderEllipsoid);

            placeholderEllipticalCylinder = BreastPhantom.addPoseToParameters( ...
                struct('a_mm', 1, 'b_mm', 1, 'length_mm', 1), ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.fatInner = AnalyticalEllipticalCylinder3D([], placeholderEllipticalCylinder);
            obj.fatOuter = AnalyticalEllipticalCylinder3D([], placeholderEllipticalCylinder);

            obj.fatComposite = CompositeAnalyticalShape3D(obj.fatOuter, obj.fatInner, params.fatIntensity);
            obj.tissueComposite = CompositeAnalyticalShape3D(obj.fatInner, [obj.beatingHeart, obj.leftLung, obj.rightLung], ...
               params.tissueIntensity);
            
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3)));
            obj.thorax = MultipleMaterialPhantom([obj.beatingHeart, obj.leftLung, obj.rightLung, obj.fatComposite, obj.tissueComposite], ...
                thoraxPose);

            breastParams = BreastPhantom.addPoseToParameters( ...
                struct('radius_mm', 1, 'length_mm', 1), ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.breastRight = AnalyticalCylinder3D([], breastParams);
            obj.breastLeft = AnalyticalCylinder3D([], breastParams);

            obj.leftAndRightBreastTissue = CompositeAnalyticalShape3D([obj.breastRight, obj.breastLeft], ...
                AnalyticalShape3D.empty(1,0), params.breastIntensity, breastParams);
            
            lesionParams = BreastPhantom.addPoseToParameters( ...
                struct('a_mm', 1, 'b_mm', 1, 'c_mm', 1), ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.lesionRight = AnalyticalEllipsoid3D(params.lesionIntensityFunction(0), ...
                lesionParams);

            obj.setBreastGeometryFromParams(params);
            obj.setShapes([obj.thorax obj.leftAndRightBreastTissue obj.lesionRight]);

            obj.updateShapesForTime(0);
        end

        function updateShapesForTime(obj, t_s)
            % updateShapesForTime  Update all sub-shape parameters at time t_s.
            %   updateShapesForTime(obj, t_s) recomputes geometry for the
            %   specified time(s) in seconds.
            validateattributes(t_s, {'numeric'}, {'real', 'finite', 'vector'}, ...
                mfilename, 't_s');

            params = obj.phantomParams;
            
            %% Heart
            centered = [0, 0, 0];
            notRotated = [0, 0, 0];
            heartParams = cardiac_ellipsoid_waveform(t_s, params.cardiacOpts);
            heartParams = BreastPhantom.addPoseToParameters(heartParams, ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.beatingHeart.setShapeParameters(heartParams);

            %% Lungs
            [lungRadius_mm, lungHeight_mm] = lung_ellipsoid_waveform(t_s, params.pulmonaryOpts);

            lungSeparation_mm = heartParams.a_mm + params.heartWallThickness_mm;
            lungPosition_mm = lungRadius_mm + lungSeparation_mm;
            lungShapeParams = struct('a_mm', lungRadius_mm, 'b_mm', lungRadius_mm, 'c_mm', lungHeight_mm);
            rightLungParams = BreastPhantom.addPoseToParameters(lungShapeParams, ...
                lungPosition_mm, zeros(size(lungPosition_mm)), zeros(size(lungPosition_mm)), ...
                notRotated(1),notRotated(2),notRotated(3));
            
            leftLungParams = BreastPhantom.addPoseToParameters(lungShapeParams, ...
                -lungPosition_mm, zeros(size(lungPosition_mm)), zeros(size(lungPosition_mm)), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.rightLung.setShapeParameters(rightLungParams);
            obj.leftLung.setShapeParameters(leftLungParams);

            %% Thorax and fat
            heart_ap_mm = heartParams.b_mm;
            chest_ap_inner_mm = max(heart_ap_mm, lungRadius_mm) + params.tissueGap_lr_mm;
            chest_lr_inner_mm = 2 * lungRadius_mm + lungSeparation_mm + params.tissueGap_lr_mm;

            fatInnerParams = struct('a_mm', chest_lr_inner_mm, ...
                'b_mm', chest_ap_inner_mm, 'length_mm', params.phantomDepth_mm);
            fatInnerParams = BreastPhantom.addPoseToParameters(fatInnerParams, ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.fatInner.setShapeParameters(fatInnerParams);

            fatThickness_mm = 10;
            chest_ap_outer_mm = chest_ap_inner_mm + fatThickness_mm;
            chest_lr_outer_mm = chest_lr_inner_mm + fatThickness_mm;
            fatOuterParams = struct('a_mm', chest_lr_outer_mm, ...
                'b_mm', chest_ap_outer_mm, 'length_mm', params.phantomDepth_mm);
            fatOuterParams = BreastPhantom.addPoseToParameters(fatOuterParams, ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.fatOuter.setShapeParameters(fatOuterParams);

            thoraxCenter = [zeros(size(chest_ap_outer_mm)), -chest_ap_outer_mm, zeros(size(chest_ap_outer_mm))];
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(...
                zeros(size(chest_ap_outer_mm)), -chest_ap_outer_mm, zeros(size(chest_ap_outer_mm)), ...
                notRotated(1),notRotated(2),notRotated(3)));
            obj.thorax.setShapeParameters(thoraxPose);

            %% Breasts
            lesionIntensity = params.lesionIntensityFunction(t_s);
            BreastPhantom.validateLesionIntensity(lesionIntensity, t_s);
            obj.lesionRight.setIntensity(lesionIntensity);
            
        end

        function S = kspaceAtTime(obj, kx, ky, kz, t_s, maxChunkSize)
            numSamples = numel(t_s);
            if ~isequal(size(kx), size(ky), size(kz), size(t_s))
                error('BreastPhantom:KspaceAtTimeSizeMismatch', ...
                    'kx, ky, kz, and t_s must have identical sizes.');
            end
            S = zeros(size(kx));

            numChunks = ceil(numSamples / maxChunkSize);
            for idxStart = 1:maxChunkSize:numSamples
                chunkIndex = ceil(idxStart / maxChunkSize);
                idxEnd = min(idxStart + maxChunkSize - 1, numSamples);
                idx = idxStart:idxEnd;
                percentComplete = 100 * chunkIndex / numChunks;
                if(numChunks>1)
                    fprintf(['kspaceAtTime: chunk %d of %d, %.1f%% complete.\n'], ...
                        chunkIndex, numChunks, percentComplete);
                end
                obj.updateShapesForTime(t_s(idx));
                S(idx) = obj.kspace(kx(idx), ky(idx), kz(idx));
            end
        end
    end

    methods (Access = private)
        function setBreastGeometryFromParams(obj, params)
            xRightBreast = params.breast_radius_mm + 0.5 * params.breast_gap_mm;
            breastParams = struct('radius_mm', params.breast_radius_mm, 'length_mm', params.breast_depth_mm);
            breastParamsRight = BreastPhantom.addPoseToParameters(breastParams, ...
                xRightBreast, 0, 0, ...
                0, 90, 90);
            breastParamsLeft = BreastPhantom.addPoseToParameters(breastParams, ...
                -xRightBreast, 0, 0, ...
                0, 90, 90);
            obj.breastRight.setShapeParameters(breastParamsRight);
            obj.breastLeft.setShapeParameters(breastParamsLeft);

            breastParamsBoth = BreastPhantom.addPoseToParameters(breastParams, ...
                0, 0.5 * params.breast_depth_mm, 0, ...
                0, 0, 0);
            obj.leftAndRightBreastTissue.setShapeParameters(breastParamsBoth);
            lesionParams = struct('a_mm', params.lesionRadius_mm, ...
                'b_mm', params.lesionRadius_mm, ...
                'c_mm', params.lesionRadius_mm);
            lesionParams = BreastPhantom.addPoseToParameters(lesionParams, ...
                xRightBreast, 0.5 * params.breast_depth_mm, 0, ...
                0, 0, 0);
            obj.lesionRight.setShapeParameters(lesionParams);
        end
    end

    methods (Access = private, Static)
        function params = addPoseToParameters(params, x_mm, y_mm, z_mm, roll, pitch, yaw)
            params.pose = BreastPhantom.createPoseStruct(x_mm, y_mm, z_mm, roll, pitch, yaw);
        end

        function pose = createPoseStruct(x_mm, y_mm, z_mm, roll, pitch, yaw)
            pose = struct('center', struct('x_mm', x_mm, ...
                'y_mm', y_mm, ...
                'z_mm', z_mm), ...
                'roll_deg', roll, ...
                'pitch_deg', pitch, ...
                'yaw_deg', yaw);
        end

        function params = normalizeLesionParams(params)
            if ~isfield(params, 'lesionRadius_mm') || isempty(params.lesionRadius_mm)
                params.lesionRadius_mm = 10;
            end
            validateattributes(params.lesionRadius_mm, {'numeric'}, ...
                {'real', 'finite', 'scalar', 'positive'}, mfilename, 'lesionRadius_mm');

            if ~isfield(params, 'lesionIntensityFunction') || isempty(params.lesionIntensityFunction)
                params.lesionIntensityFunction = @(t_s) min(2, max(0, 2 .* t_s ./ 100));
            end
            if ~isa(params.lesionIntensityFunction, 'function_handle')
                error('BreastPhantom:InvalidLesionIntensityFunction', ...
                    'lesionIntensityFunction must be a function handle.');
            end
        end

        function validateLesionIntensity(intensity, t_s)
            if ~isnumeric(intensity)
                error('BreastPhantom:InvalidLesionIntensity', ...
                    'Lesion intensity must be numeric.');
            end
            if ~(isscalar(intensity) || isequal(size(intensity), size(t_s)))
                error('BreastPhantom:InvalidLesionIntensitySize', ...
                    'Lesion intensity must be scalar or the same size as t_s.');
            end
            validateattributes(intensity, {'numeric'}, {'real', 'finite'}, ...
                mfilename, 'lesionIntensity');
        end
    end
end
