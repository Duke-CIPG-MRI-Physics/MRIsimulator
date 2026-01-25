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
                struct('a_mm', 1, 'b_mm', 1, 'c_mm', 1), centered, notRotated);
            obj.beatingHeart = AnalyticalEllipsoid3D(params.heartIntensity, placeholderEllipsoid);
            obj.rightLung = AnalyticalEllipsoid3D(params.lungIntensity, placeholderEllipsoid);
            obj.leftLung = AnalyticalEllipsoid3D(params.lungIntensity, placeholderEllipsoid);

            placeholderEllipticalCylinder = BreastPhantom.addPoseToParameters( ...
                struct('a_mm', 1, 'b_mm', 1, 'length_mm', 1), centered, notRotated);
            obj.fatInner = AnalyticalEllipticalCylinder3D([], placeholderEllipticalCylinder);
            obj.fatOuter = AnalyticalEllipticalCylinder3D([], placeholderEllipticalCylinder);

            obj.fatComposite = CompositeAnalyticalShape3D(obj.fatOuter, obj.fatInner, params.fatIntensity);
            obj.tissueComposite = CompositeAnalyticalShape3D(obj.fatInner, [obj.beatingHeart, obj.leftLung, obj.rightLung], ...
               params.tissueIntensity);
            
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(centered, notRotated));
            obj.thorax = MultipleMaterialPhantom([obj.beatingHeart, obj.leftLung, obj.rightLung, obj.fatComposite, obj.tissueComposite], ...
                thoraxPose);

            placeholderCylinder = BreastPhantom.addPoseToParameters( ...
                struct('radius_mm', 1, 'length_mm', 1), centered, notRotated);
            obj.breastRight = AnalyticalCylinder3D([], placeholderCylinder);
            obj.breastLeft = AnalyticalCylinder3D([], placeholderCylinder);

            obj.enhancedCylinder = AnalyticalCylinder3D(params.enhancedIntensity, placeholderCylinder);
            obj.unenhancedCylinder = AnalyticalCylinder3D(params.unenhancedIntensity, placeholderCylinder);

            obj.leftAndRightBreastTissue = CompositeAnalyticalShape3D([obj.breastRight, obj.breastLeft], ...
                [obj.unenhancedCylinder, obj.enhancedCylinder], params.breastIntensity, placeholderCylinder);

            obj.enhancingVessel = MultipleMaterialPhantom([obj.enhancedCylinder, obj.unenhancedCylinder], placeholderCylinder);

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
            heartParams = BreastPhantom.addPoseToParameters(heartParams, centered, notRotated);
            obj.beatingHeart.setShapeParameters(heartParams);

            %% Lungs
            [lungRadius_mm, lungHeight_mm] = lung_ellipsoid_waveform(t_s, params.pulmonaryOpts);
            lungRadius_mm = lungRadius_mm(:);
            lungHeight_mm = lungHeight_mm(:);

            lungSeparation_mm = heartParams.a_mm + params.heartWallThickness_mm;
            lungPosition_mm = lungRadius_mm + lungSeparation_mm;

            lungShapeParams = struct('a_mm', lungRadius_mm, 'b_mm', lungRadius_mm, 'c_mm', lungHeight_mm);
            rightLungCenter = [lungPosition_mm, zeros(size(lungPosition_mm)), zeros(size(lungPosition_mm))];
            leftLungCenter = [-lungPosition_mm, zeros(size(lungPosition_mm)), zeros(size(lungPosition_mm))];

            rightLungParams = BreastPhantom.addPoseToParameters(lungShapeParams, rightLungCenter, notRotated);
            leftLungParams = BreastPhantom.addPoseToParameters(lungShapeParams, leftLungCenter, notRotated);
            obj.rightLung.setShapeParameters(rightLungParams);
            obj.leftLung.setShapeParameters(leftLungParams);

            %% Thorax and fat
            heart_ap_mm = heartParams.b_mm;
            chest_ap_inner_mm = max(heart_ap_mm, lungRadius_mm) + params.tissueGap_lr_mm;
            chest_lr_inner_mm = 2 * lungRadius_mm + lungSeparation_mm + params.tissueGap_lr_mm;

            fatInnerParams = struct('a_mm', chest_lr_inner_mm, ...
                'b_mm', chest_ap_inner_mm, 'length_mm', params.phantomDepth_mm);
            fatInnerParams = BreastPhantom.addPoseToParameters(fatInnerParams, centered, notRotated);
            obj.fatInner.setShapeParameters(fatInnerParams);

            fatThickness_mm = 10;
            chest_ap_outer_mm = chest_ap_inner_mm + fatThickness_mm;
            chest_lr_outer_mm = chest_lr_inner_mm + fatThickness_mm;
            fatOuterParams = struct('a_mm', chest_lr_outer_mm, ...
                'b_mm', chest_ap_outer_mm, 'length_mm', params.phantomDepth_mm);
            fatOuterParams = BreastPhantom.addPoseToParameters(fatOuterParams, centered, notRotated);
            obj.fatOuter.setShapeParameters(fatOuterParams);

            thoraxCenter = [zeros(size(chest_ap_outer_mm)), -chest_ap_outer_mm, zeros(size(chest_ap_outer_mm))];
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(thoraxCenter, notRotated));
            obj.thorax.setShapeParameters(thoraxPose);

            %% Breasts
            [rightBreastCenter, leftBreastCenter, breastCenter] = obj.breastCenters(params);
            breastParams = struct('radius_mm', params.breast_radius_mm, 'length_mm', params.breast_depth_mm);
            breastParamsRight = BreastPhantom.addPoseToParameters(breastParams, rightBreastCenter, [0, 90, 90]);
            breastParamsLeft = BreastPhantom.addPoseToParameters(breastParams, leftBreastCenter, [0, 90, 90]);
            obj.breastRight.setShapeParameters(breastParamsRight);
            obj.breastLeft.setShapeParameters(breastParamsLeft);

            breastParamsBoth = BreastPhantom.addPoseToParameters(breastParams, breastCenter, notRotated);
            obj.leftAndRightBreastTissue.setShapeParameters(breastParamsBoth);
            obj.enhancingVessel.setShapeParameters(breastParamsBoth);

            %% Vessel enhancement
            enhancedParams = obj.vesselParameters('enhanced', t_s, params, rightBreastCenter);
            unenhancedParams = obj.vesselParameters('unenhanced', t_s, params, rightBreastCenter);
            obj.enhancedCylinder.setShapeParameters(enhancedParams);
            obj.unenhancedCylinder.setShapeParameters(unenhancedParams);
        end

        function S = kspaceAtTime(obj, kx, ky, kz, t_s)
            maxChunkSize = 15625;
            numSamples = numel(t_s);
            if ~isequal(size(kx), size(ky), size(kz), size(t_s))
                error('BreastPhantom:KspaceAtTimeSizeMismatch', ...
                    'kx, ky, kz, and t_s must be the same size.');
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

    methods (Access = private, Static)
        function params = addPoseToParameters(params, centerVec, rollPitchYaw)
            params.pose = BreastPhantom.createPoseStruct(centerVec, rollPitchYaw);
        end

        function pose = createPoseStruct(centerVec, rollPitchYaw)
            if size(centerVec, 1) == 3 && size(centerVec, 2) ~= 3
                x_mm = centerVec(1, :);
                y_mm = centerVec(2, :);
                z_mm = centerVec(3, :);
            elseif size(centerVec, 2) == 3
                x_mm = centerVec(:, 1);
                y_mm = centerVec(:, 2);
                z_mm = centerVec(:, 3);
            else
                error('BreastPhantom:InvalidCenterVec', ...
                    'centerVec must be 3xN (axis rows) or Nx3 (axis columns).');
            end

            pose = struct('center', struct('x_mm', x_mm, ...
                'y_mm', y_mm, ...
                'z_mm', z_mm), ...
                'roll_deg', rollPitchYaw(1), ...
                'pitch_deg', rollPitchYaw(2), ...
                'yaw_deg', rollPitchYaw(3));
        end
    end

    methods (Access = private)
        function vesselParams = vesselParameters(~, segment, t_s, params, rightBreastCenter)
            elapsedTime_s = max(t_s - params.startInjectionTime_s, 0);
            enhancedLength_mm = min(params.breastVesselVelocity_cm_s * 10 .* elapsedTime_s, ...
                params.totalVesselLength_mm);
            unenhancedLength_mm = max(params.totalVesselLength_mm - enhancedLength_mm, 0);

            switch lower(segment)
                case 'enhanced'
                    length_mm = enhancedLength_mm(:);
                    center_mm = repmat(rightBreastCenter, numel(length_mm), 1);
                    center_mm(:, 3) = center_mm(:, 3) - 0.5 .* unenhancedLength_mm(:);
                case 'unenhanced'
                    length_mm = unenhancedLength_mm(:);
                    center_mm = repmat(rightBreastCenter, numel(length_mm), 1);
                    center_mm(:, 3) = center_mm(:, 3) + 0.5 .* enhancedLength_mm(:);
                otherwise
                    error('BreastPhantom:InvalidVesselSegment', ...
                        'Segment must be ''enhanced'' or ''unenhanced''.');
            end

            vesselParams = struct('radius_mm', params.vesselRadius_mm, 'length_mm', length_mm);
            vesselParams.pose = BreastPhantom.createPoseStruct(center_mm, params.breastRollPitchYaw);
        end

        function [rightBreastCenter, leftBreastCenter, breastCenter] = breastCenters(~, params)
            rightBreastCenter = [params.breast_radius_mm + 0.5 * params.breast_gap_mm, 0, 0];
            leftBreastCenter = [-rightBreastCenter(1), rightBreastCenter(2), rightBreastCenter(3)];
            breastCenter = [0, 0.5 * params.breast_depth_mm, 0];
        end
    end

end
