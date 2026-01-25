classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
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
        leftAndRightBreastTissue
        enhancedCylinder
        unenhancedCylinder
        enhancingVessel
    end

    methods
        function obj = BreastPhantom(params)
            arguments
                params struct = createBreastPhantomParams();
            end

            obj@MultipleMaterialPhantom();
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

            placeholderCylinder = BreastPhantom.addPoseToParameters( ...
                struct('radius_mm', 1, 'length_mm', 1), ...
                centered(1),centered(2),centered(3), ...
                notRotated(1),notRotated(2),notRotated(3));
            obj.breastRight = AnalyticalCylinder3D([], placeholderCylinder);
            obj.breastLeft = AnalyticalCylinder3D([], placeholderCylinder);

            % Use vessel intensity deltas so we can avoid subtracting and re-adding
            % identical cylinder k-space when combining with breast tissue.
            enhancedDelta = params.enhancedIntensity - params.breastIntensity;
            unenhancedDelta = params.unenhancedIntensity - params.breastIntensity;
            obj.enhancedCylinder = AnalyticalCylinder3D(enhancedDelta, placeholderCylinder);
            obj.unenhancedCylinder = AnalyticalCylinder3D(unenhancedDelta, placeholderCylinder);

            obj.leftAndRightBreastTissue = CompositeAnalyticalShape3D([obj.breastRight, obj.breastLeft], ...
                AnalyticalShape3D.empty(1,0), params.breastIntensity, placeholderCylinder);

            obj.enhancingVessel = MultipleMaterialPhantom([obj.enhancedCylinder, obj.unenhancedCylinder], ...
                placeholderCylinder);

            % enhancingVessel  = CompositeAnalyticalShape3D([unenhancedCylinder, enhancedCylinder], [], ...
            %     breastIntensity , breastParamsBoth);

            obj.setShapes([obj.thorax obj.leftAndRightBreastTissue obj.enhancingVessel]);

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
                notRotated(1),notRotated(2),notRotated(3));
            obj.leftAndRightBreastTissue.setShapeParameters(breastParamsBoth);
            obj.enhancingVessel.setShapeParameters(breastParamsBoth);


    
            timePostInj_s = max(t_s - params.startInjectionTime_s, 0);
            enhancedLength_mm = min(params.breastVesselVelocity_cm_s * 10 .* timePostInj_s, ...
                params.totalVesselLength_mm);
            unenhancedLength_mm = max(params.totalVesselLength_mm - enhancedLength_mm, 0);

            rightBreastCenter_z = 0;

            % case 'enhanced'
            center_enhanced_mm = rightBreastCenter_z - 0.5 .* unenhancedLength_mm;
            % case 'unenhanced'
            center_unenhanced_mm = rightBreastCenter_z + 0.5 .* enhancedLength_mm;

            enhancedVesselParams = struct('radius_mm', params.vesselRadius_mm, ...
                'length_mm', enhancedLength_mm);
            enhancedVesselParams.pose = BreastPhantom.createPoseStruct(...
                0, 0, center_enhanced_mm, ...
                params.breastRollPitchYaw(1), params.breastRollPitchYaw(2), params.breastRollPitchYaw(3));

            unenhancedVesselParams = struct('radius_mm', params.vesselRadius_mm, ...
                'length_mm', unenhancedLength_mm);
            unenhancedVesselParams.pose = BreastPhantom.createPoseStruct(...
                0, 0, center_unenhanced_mm, ...
                params.breastRollPitchYaw(1), params.breastRollPitchYaw(2), params.breastRollPitchYaw(3));

            %% Vessel enhancement
            obj.enhancedCylinder.setShapeParameters(enhancedVesselParams);
            obj.unenhancedCylinder.setShapeParameters(unenhancedVesselParams);
        end

        function S = kspaceAtTime(obj, kx, ky, kz, t_s)
            maxChunkSize = 250000;
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
                fprintf(['kspaceAtTime: chunk %d of %d, %.1f%% complete.\n'], ...
                    chunkIndex, numChunks, percentComplete);
                obj.updateShapesForTime(t_s(idx));
                S(idx) = obj.kspace(kx(idx), ky(idx), kz(idx));
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
    end
end
