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
        enhancedCylinders
        unenhancedCylinders
        enhancingVessel
    end

    methods
        function obj = BreastPhantom(params)
            arguments
                params struct = createBreastPhantomParams();
            end

            obj@MultipleMaterialPhantom();
            params = BreastPhantom.normalizeVesselParams(params);
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
            numSegments = numel(params.vesselSegments);
            obj.enhancedCylinders = AnalyticalCylinder3D.empty(1, 0);
            obj.unenhancedCylinders = AnalyticalCylinder3D.empty(1, 0);
            for idx = 1:numSegments
                obj.enhancedCylinders(idx) = AnalyticalCylinder3D(enhancedDelta, placeholderCylinder);
                obj.unenhancedCylinders(idx) = AnalyticalCylinder3D(unenhancedDelta, placeholderCylinder);
            end

            obj.leftAndRightBreastTissue = CompositeAnalyticalShape3D([obj.breastRight, obj.breastLeft], ...
                AnalyticalShape3D.empty(1,0), params.breastIntensity, placeholderCylinder);

            obj.enhancingVessel = MultipleMaterialPhantom([obj.enhancedCylinders, obj.unenhancedCylinders], ...
                placeholderCylinder);

            % enhancingVessel  = CompositeAnalyticalShape3D([unenhancedCylinder, enhancedCylinder], [], ...
            %     breastIntensity , breastParamsBoth);

            obj.setBreastGeometryFromParams(params);
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
            %% Vessel enhancement
            timePostInj_s = max(t_s - params.startInjectionTime_s, 0);
            flowLength_mm = params.breastVesselVelocity_cm_s * 10 .* timePostInj_s;
            cumulativeLength_mm = 0;

            for idx = 1:numel(params.vesselSegments)
                segment = params.vesselSegments(idx);
                segmentLength_mm = segment.length_mm;
                enhancedLength_mm = min(max(flowLength_mm - cumulativeLength_mm, 0), segmentLength_mm);
                unenhancedLength_mm = segmentLength_mm - enhancedLength_mm;
                cumulativeLength_mm = cumulativeLength_mm + segmentLength_mm;

                [segmentCenter_mm, axisUnit] = BreastPhantom.segmentCenterAndAxis(segment.pose);
                segmentStart_mm = segmentCenter_mm - axisUnit .* (0.5 * segmentLength_mm);
                enhancedCenter_mm = segmentStart_mm + axisUnit .* (0.5 * enhancedLength_mm);
                unenhancedCenter_mm = segmentStart_mm + axisUnit .* (enhancedLength_mm + 0.5 * unenhancedLength_mm);

                enhancedVesselParams = struct('radius_mm', segment.radius_mm, ...
                    'length_mm', enhancedLength_mm);
                enhancedVesselParams.pose = BreastPhantom.createPoseStruct(...
                    enhancedCenter_mm(1, :), enhancedCenter_mm(2, :), enhancedCenter_mm(3, :), ...
                    segment.pose.roll_deg, segment.pose.pitch_deg, segment.pose.yaw_deg);

                unenhancedVesselParams = struct('radius_mm', segment.radius_mm, ...
                    'length_mm', unenhancedLength_mm);
                unenhancedVesselParams.pose = BreastPhantom.createPoseStruct(...
                    unenhancedCenter_mm(1, :), unenhancedCenter_mm(2, :), unenhancedCenter_mm(3, :), ...
                    segment.pose.roll_deg, segment.pose.pitch_deg, segment.pose.yaw_deg);

                obj.enhancedCylinders(idx).setShapeParameters(enhancedVesselParams);
                obj.unenhancedCylinders(idx).setShapeParameters(unenhancedVesselParams);
            end
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
            obj.enhancingVessel.setShapeParameters(breastParamsBoth);
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

        function params = normalizeVesselParams(params)
            if ~isfield(params, 'vesselSegments') || isempty(params.vesselSegments)
                params.vesselSegments = BreastPhantom.buildDefaultVesselSegments(params);
            end
            params.vesselSegments = BreastPhantom.ensureVesselSegmentDefaults( ...
                params.vesselSegments, params);
        end

        function segments = buildDefaultVesselSegments(params)
            segmentCount = params.vesselSegmentCount;
            segmentLength_mm = params.totalVesselLength_mm / segmentCount;
            segmentLengths_mm = repmat(segmentLength_mm, 1, segmentCount);
            segments = createBreastVesselSegments(segmentLengths_mm, params.vesselRadius_mm, ...
                [0 0 0], params.breastRollPitchYaw);
        end

        function segments = ensureVesselSegmentDefaults(segments, params)
            defaultPose = BreastPhantom.createPoseStruct(0, 0, 0, ...
                params.breastRollPitchYaw(1), params.breastRollPitchYaw(2), params.breastRollPitchYaw(3));

            for idx = 1:numel(segments)
                if ~isfield(segments(idx), 'length_mm') || isempty(segments(idx).length_mm)
                    error('BreastPhantom:MissingVesselSegmentLength', ...
                        'Each vessel segment must define length_mm in mm.');
                end
                validateattributes(segments(idx).length_mm, {'numeric'}, ...
                    {'real', 'finite', 'scalar', 'positive'}, mfilename, ...
                    sprintf('vesselSegments(%d).length_mm', idx));

                if ~isfield(segments(idx), 'radius_mm') || isempty(segments(idx).radius_mm)
                    segments(idx).radius_mm = params.vesselRadius_mm;
                end
                validateattributes(segments(idx).radius_mm, {'numeric'}, ...
                    {'real', 'finite', 'scalar', 'positive'}, mfilename, ...
                    sprintf('vesselSegments(%d).radius_mm', idx));

                if ~isfield(segments(idx), 'pose') || isempty(segments(idx).pose)
                    segments(idx).pose = defaultPose;
                else
                    segments(idx).pose = BreastPhantom.mergePoseDefaults( ...
                        segments(idx).pose, defaultPose);
                end
            end
        end

        function pose = mergePoseDefaults(pose, defaultPose)
            if ~isfield(pose, 'center') || ~isstruct(pose.center)
                pose.center = defaultPose.center;
            else
                if ~isfield(pose.center, 'x_mm')
                    pose.center.x_mm = defaultPose.center.x_mm;
                end
                if ~isfield(pose.center, 'y_mm')
                    pose.center.y_mm = defaultPose.center.y_mm;
                end
                if ~isfield(pose.center, 'z_mm')
                    pose.center.z_mm = defaultPose.center.z_mm;
                end
            end

            if ~isfield(pose, 'roll_deg')
                pose.roll_deg = defaultPose.roll_deg;
            end
            if ~isfield(pose, 'pitch_deg')
                pose.pitch_deg = defaultPose.pitch_deg;
            end
            if ~isfield(pose, 'yaw_deg')
                pose.yaw_deg = defaultPose.yaw_deg;
            end
        end

        function [center_mm, axisUnit] = segmentCenterAndAxis(pose)
            center_mm = [pose.center.x_mm; pose.center.y_mm; pose.center.z_mm];
            axisUnit = BreastPhantom.axisUnitFromPose(pose);
        end

        function axisUnit = axisUnitFromPose(pose)
            roll_rad = deg2rad(pose.roll_deg);
            pitch_rad = deg2rad(pose.pitch_deg);
            yaw_rad = deg2rad(pose.yaw_deg);

            cr = cos(roll_rad);
            sr = sin(roll_rad);
            cp = cos(pitch_rad);
            sp = sin(pitch_rad);
            cy = cos(yaw_rad);
            sy = sin(yaw_rad);

            axisUnit = [cy .* sp .* cr + sy .* sr; ...
                sy .* sp .* cr - cy .* sr; ...
                cp .* cr];
        end
    end
end
