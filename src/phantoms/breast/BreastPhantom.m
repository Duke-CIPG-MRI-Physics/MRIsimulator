classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and
    %   enhancing lesions in the breast tissue. Lesion center coordinates
    %   are specified relative to the center of the right breast and are
    %   grouped with that breast in a nested MultipleMaterialPhantom.
    %   Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup. Heart geometry leverages
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
        rightBreastPhantom
        lesions
    end

    methods
        function obj = BreastPhantom(params)
            arguments
                params struct = createBreastPhantomParams();
            end

            obj@MultipleMaterialPhantom();
            params = BreastPhantom.normalizeLesionParams(params);
            obj.phantomParams = params;

            centered = [0, 0, 0];
            notRotated = [0, 0, 0];

            placeholderEllipsoid = BreastPhantom.addPoseToParameters( ...
                struct('a_mm', 1, 'b_mm', 1, 'c_mm', 1), ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.beatingHeart = AnalyticalEllipsoid3D(params.heartIntensity, placeholderEllipsoid);
            obj.rightLung = AnalyticalEllipsoid3D(params.lungIntensity, placeholderEllipsoid);
            obj.leftLung = AnalyticalEllipsoid3D(params.lungIntensity, placeholderEllipsoid);

            placeholderEllipticalCylinder = BreastPhantom.addPoseToParameters( ...
                struct('a_mm', 1, 'b_mm', 1, 'length_mm', 1), ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.fatInner = AnalyticalEllipticalCylinder3D([], placeholderEllipticalCylinder);
            obj.fatOuter = AnalyticalEllipticalCylinder3D([], placeholderEllipticalCylinder);

            obj.fatComposite = CompositeAnalyticalShape3D(obj.fatOuter, obj.fatInner, params.fatIntensity);
            obj.tissueComposite = CompositeAnalyticalShape3D(obj.fatInner, [obj.beatingHeart, obj.leftLung, obj.rightLung], ...
                params.tissueIntensity);

            thoraxPose = struct('pose', BreastPhantom.createPoseStruct( ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3)));
            obj.thorax = MultipleMaterialPhantom( ...
                [obj.beatingHeart, obj.leftLung, obj.rightLung, obj.fatComposite, obj.tissueComposite], ...
                thoraxPose);

            breastParams = BreastPhantom.addPoseToParameters( ...
                struct('radius_mm', 1, 'length_mm', 1), ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.breastRight = AnalyticalCylinder3D(params.breastIntensity, breastParams);
            obj.breastLeft = AnalyticalCylinder3D(params.breastIntensity, breastParams);

            obj.lesions = AnalyticalShape3D.empty(1, 0);
            for lesionIdx = 1:numel(params.lesions)
                lesionDef = params.lesions(lesionIdx);
                lesionParams = BreastPhantom.createLesionShapeParameters(lesionDef);
                obj.lesions(lesionIdx) = AnalyticalEllipsoid3D(lesionDef.intensityFunction(0), lesionParams);
            end

            rightBreastPose = struct('pose', BreastPhantom.createPoseStruct( ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3)));
            obj.rightBreastPhantom = MultipleMaterialPhantom([obj.breastRight, obj.lesions], rightBreastPose);

            obj.setBreastGeometryFromParams(params);
            obj.setShapes([obj.thorax, obj.breastLeft, obj.rightBreastPhantom]);
        end

        function updateShapesForTime(obj, t_s)
            validateattributes(t_s, {'numeric'}, {'real', 'finite', 'vector'}, ...
                mfilename, 't_s');

            params = obj.phantomParams;
            centered = [0, 0, 0];
            notRotated = [0, 0, 0];

            heartParams = cardiac_ellipsoid_waveform(t_s, params.cardiacOpts);
            heartParams = BreastPhantom.addPoseToParameters(heartParams, ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.beatingHeart.setShapeParameters(heartParams);

            [lungRadius_mm, lungHeight_mm] = lung_ellipsoid_waveform(t_s, params.pulmonaryOpts);

            lungSeparation_mm = heartParams.a_mm + params.heartWallThickness_mm;
            lungPosition_mm = lungRadius_mm + lungSeparation_mm;
            lungShapeParams = struct('a_mm', lungRadius_mm, 'b_mm', lungRadius_mm, 'c_mm', lungHeight_mm);
            rightLungParams = BreastPhantom.addPoseToParameters(lungShapeParams, ...
                lungPosition_mm, zeros(size(lungPosition_mm)), zeros(size(lungPosition_mm)), ...
                notRotated(1), notRotated(2), notRotated(3));
            leftLungParams = BreastPhantom.addPoseToParameters(lungShapeParams, ...
                -lungPosition_mm, zeros(size(lungPosition_mm)), zeros(size(lungPosition_mm)), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.rightLung.setShapeParameters(rightLungParams);
            obj.leftLung.setShapeParameters(leftLungParams);

            heart_ap_mm = heartParams.b_mm;
            chest_ap_inner_mm = max(heart_ap_mm, lungRadius_mm) + params.tissueGap_lr_mm;
            chest_lr_inner_mm = 2 * lungRadius_mm + lungSeparation_mm + params.tissueGap_lr_mm;

            fatInnerParams = struct('a_mm', chest_lr_inner_mm, ...
                'b_mm', chest_ap_inner_mm, 'length_mm', params.phantomDepth_mm);
            fatInnerParams = BreastPhantom.addPoseToParameters(fatInnerParams, ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.fatInner.setShapeParameters(fatInnerParams);

            fatThickness_mm = 10;
            chest_ap_outer_mm = chest_ap_inner_mm + fatThickness_mm;
            chest_lr_outer_mm = chest_lr_inner_mm + fatThickness_mm;
            fatOuterParams = struct('a_mm', chest_lr_outer_mm, ...
                'b_mm', chest_ap_outer_mm, 'length_mm', params.phantomDepth_mm);
            fatOuterParams = BreastPhantom.addPoseToParameters(fatOuterParams, ...
                centered(1), centered(2), centered(3), ...
                notRotated(1), notRotated(2), notRotated(3));
            obj.fatOuter.setShapeParameters(fatOuterParams);

            thoraxPose = struct('pose', BreastPhantom.createPoseStruct( ...
                zeros(size(chest_ap_outer_mm)), -chest_ap_outer_mm, zeros(size(chest_ap_outer_mm)), ...
                notRotated(1), notRotated(2), notRotated(3)));
            obj.thorax.setShapeParameters(thoraxPose);

            for lesionIdx = 1:numel(obj.lesions)
                lesionIntensity = params.lesions(lesionIdx).intensityFunction(t_s);
                BreastPhantom.validateLesionIntensity(lesionIntensity, t_s);
                obj.lesions(lesionIdx).setIntensity(@() lesionIntensity);
            end
        end

        function S = kspaceAtTime(obj, kx, ky, kz, t_s, maxChunkSize)
            numSamples = numel(t_s);
            if ~isequal(size(kx), size(ky), size(kz), size(t_s))
                error('BreastPhantom:KspaceAtTimeSizeMismatch', ...
                    'kx, ky, kz, and t_s must have identical sizes.');
            end
            S = zeros(size(kx));

            for idxStart = 1:maxChunkSize:numSamples
                idxEnd = min(idxStart + maxChunkSize - 1, numSamples);
                idx = idxStart:idxEnd;
                obj.updateShapesForTime(t_s(idx));
                S(idx) = obj.kspace(kx(idx), ky(idx), kz(idx));
            end
        end
    end

    methods (Access = private)
        function setBreastGeometryFromParams(obj, params)
            xRightBreast = params.breast_radius_mm + 0.5 * params.breast_gap_mm;
            breastParams = struct('radius_mm', params.breast_radius_mm, 'length_mm', params.breast_depth_mm);
            breastParamsLeft = BreastPhantom.addPoseToParameters(breastParams, ...
                -xRightBreast, 0.5 * params.breast_depth_mm, 0, ...
                0, 90, 90);
            obj.breastLeft.setShapeParameters(breastParamsLeft);

            breastParamsRightLocal = BreastPhantom.addPoseToParameters(breastParams, ...
                0, 0, 0, ...
                0, 90, 90);
            obj.breastRight.setShapeParameters(breastParamsRightLocal);

            rightBreastPose = struct('pose', BreastPhantom.createPoseStruct( ...
                xRightBreast, 0.5 * params.breast_depth_mm, 0, ...
                0, 0, 0));
            obj.rightBreastPhantom.setShapeParameters(rightBreastPose);

            for lesionIdx = 1:numel(obj.lesions)
                obj.lesions(lesionIdx).setShapeParameters( ...
                    BreastPhantom.createLesionShapeParameters(params.lesions(lesionIdx)));
            end
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

        function lesionParams = createLesionShapeParameters(lesionDef)
            center_mm = lesionDef.center_mm;
            lesionParams = struct('a_mm', lesionDef.radius_mm, ...
                'b_mm', lesionDef.radius_mm, ...
                'c_mm', lesionDef.radius_mm);
            lesionParams = BreastPhantom.addPoseToParameters(lesionParams, ...
                center_mm(1), center_mm(2), center_mm(3), 0, 0, 0);
        end

        function params = normalizeLesionParams(params)
            if ~isfield(params, 'lesionRadius_mm') || isempty(params.lesionRadius_mm)
                params.lesionRadius_mm = 10;
            end
            validateattributes(params.lesionRadius_mm, {'numeric'}, ...
                {'real', 'finite', 'scalar', 'positive'}, mfilename, 'lesionRadius_mm');

            if ~isfield(params, 'lesionWashinType') || isempty(params.lesionWashinType)
                params.lesionWashinType = "medium";
            end
            if ~isfield(params, 'lesionWashoutType') || isempty(params.lesionWashoutType)
                params.lesionWashoutType = "plateau";
            end
            if ~isfield(params, 'lesionKineticOverrides') || isempty(params.lesionKineticOverrides)
                params.lesionKineticOverrides = struct();
            end

            if ~isfield(params, 'lesionIntensityFunction') || isempty(params.lesionIntensityFunction)
                params.lesionIntensityFunction = @(t_s) calculateLesionEnhancement( ...
                    t_s, params, params.lesionWashinType, params.lesionWashoutType, ...
                    params.lesionKineticOverrides);
            end
            if ~isa(params.lesionIntensityFunction, 'function_handle')
                error('BreastPhantom:InvalidLesionIntensityFunction', ...
                    'lesionIntensityFunction must be a function handle.');
            end

            if ~isfield(params, 'lesions') || isempty(params.lesions)
                params.lesions = struct( ...
                    'center_mm', [0, 0, 0], ...
                    'radius_mm', params.lesionRadius_mm, ...
                    'intensityFunction', params.lesionIntensityFunction);
            end

            params.lesions = BreastPhantom.normalizeLesionArray(params.lesions, params.lesionIntensityFunction);
        end

        function lesions = normalizeLesionArray(lesions, defaultIntensityFunction)
            if ~isstruct(lesions) || isempty(lesions)
                error('BreastPhantom:InvalidLesions', ...
                    'lesions must be a non-empty struct array.');
            end

            for lesionIdx = 1:numel(lesions)
                lesionDef = lesions(lesionIdx);
                if ~isfield(lesionDef, 'center_mm') || isempty(lesionDef.center_mm)
                    error('BreastPhantom:MissingLesionCenter', ...
                        'lesions(%d).center_mm must be a 1x3 vector in mm relative to the right breast center.', lesionIdx);
                end
                validateattributes(lesionDef.center_mm, {'numeric'}, ...
                    {'real', 'finite', 'numel', 3}, mfilename, sprintf('lesions(%d).center_mm', lesionIdx));

                if ~isfield(lesionDef, 'radius_mm') || isempty(lesionDef.radius_mm)
                    error('BreastPhantom:MissingLesionRadius', ...
                        'lesions(%d).radius_mm must be a positive scalar in mm.', lesionIdx);
                end
                validateattributes(lesionDef.radius_mm, {'numeric'}, ...
                    {'real', 'finite', 'scalar', 'positive'}, mfilename, sprintf('lesions(%d).radius_mm', lesionIdx));

                if ~isfield(lesionDef, 'intensityFunction') || isempty(lesionDef.intensityFunction)
                    lesions(lesionIdx).intensityFunction = defaultIntensityFunction;
                elseif ~isa(lesionDef.intensityFunction, 'function_handle')
                    error('BreastPhantom:InvalidLesionIntensityFunction', ...
                        'lesions(%d).intensityFunction must be a function handle.', lesionIdx);
                end

                lesions(lesionIdx).center_mm = reshape(lesionDef.center_mm, 1, 3);
                lesions(lesionIdx).radius_mm = lesionDef.radius_mm;
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
