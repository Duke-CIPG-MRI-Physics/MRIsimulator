classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup. Heart geometry leverages
    %   cardiac_ellipsoid_waveform, which evaluates long time vectors in
    %   ~0.001 Gb, 15,625-sample chunks to keep memory usage minimal.

    properties (Access = private)
        time_s (1,:) double {mustBeFinite}
    end

    methods
        function obj = BreastPhantom(t_s)
            arguments
                t_s (1,:) double {mustBeFinite}
            end

            obj@MultipleMaterialPhantom();

            obj.time_s = t_s;

            %% Create heart
            cardiacOpts = struct('HR_bpm', 70/10.66, ...
                'EDV_ml', 150,...
                'ESV_ml', 75,...
                'systFrac', 0.35, ...
                'q_ED', 50/27, ...
                'GLS_peak', -0.20, ...
                'GCS_peak', -0.25);
            heartIntensity = 1;
            heartWallThickness_mm = 8;
            centered = [0, 0, 0];
            notRotated = [0, 0, 0];

            cardiacParams = @()BreastPhantom.addPoseToParameters( ...
                cardiac_ellipsoid_waveform(t_s, cardiacOpts), centered, notRotated);
            beatingHeart = AnalyticalEllipsoid3D(heartIntensity, cardiacParams);

            %% Create lungs
            pulmonaryOpts = struct('f_bpm', 12/10.66, ...
                'VT_L', 0.4,...
                'Vres_L', 0.8, ...
                'Vbase_L', 1.5, ...
                'bellyFrac', 0.6, ...
                'inspFrac', 1/3, ...
                'GCS_peak', 0.4);
            lungIntensity = 0.1;

            % ROBBIE = here are the cardiac values
            % beatingHeart.getShapeParameters.a_mm -< I think a is the r/l
            % and a/p direction and b is the s/i direction 

            % TODO make this a function handle so its not directly
            % evaluated yet.
            pulmonaryOpts.lungSeparation_mm = beatingHeart.getShapeParameters.a_mm + heartWallThickness_mm; 

            lungShapeParameters = BreastPhantom.addPoseToParameters(struct(), centered, notRotated);
            breathingLung = BreathingLung(obj.time_s, pulmonaryOpts, ...
                lungIntensity, lungShapeParameters);

            %% Create thorax
            bodyShift = -80;
            tissueGap_lr_mm = 30;
            phantomDepth_mm = 300;
            fatIntensity = 2;
            tissueIntensity = 0.5;
            heart_ap_mm = beatingHeart.getShapeParameters.b_mm;
            lungRadius = breathingLung.getLeftLung().getShapeParameters().R_mm;

            chest_ap_inner_mm = max(heart_ap_mm, lungRadius) + tissueGap_lr_mm;
            chest_lr_inner_mm = 2 * lungRadius + pulmonaryOpts.lungSeparation_mm + tissueGap_lr_mm;
            
            % ROBBIE - sanity check that the chest dimmensions make sense
            % in the oversampled FOV

            fatInnerParams = struct('a_mm', chest_lr_inner_mm, ...
                'b_mm', chest_ap_inner_mm, 'length_mm', phantomDepth_mm);
            fatInnerParams = BreastPhantom.addPoseToParameters(fatInnerParams, [0, 0, 0], [0, 0, 0]);
            fat_inner = AnalyticalEllipticalCylinder3D([], fatInnerParams);

            fatThickness_mm = 10;
            chest_ap_outer_mm = chest_ap_inner_mm + fatThickness_mm;
            chest_lr_outer_mm = chest_lr_inner_mm + fatThickness_mm;
            fatOuterParams = struct('a_mm', chest_lr_outer_mm, ...
                'b_mm', chest_ap_outer_mm, 'length_mm', phantomDepth_mm);
            fatOuterParams = BreastPhantom.addPoseToParameters(fatOuterParams, [0, 0, 0], [0, 0, 0]);
            fat_outer = AnalyticalEllipticalCylinder3D([], fatOuterParams);

            fatComposite = CompositeAnalyticalShape3D(fat_outer, fat_inner, fatIntensity);
            tissueComposite = CompositeAnalyticalShape3D(fat_inner, [beatingHeart, breathingLung], ...
               tissueIntensity);
            
            thoraxCenter = [zeros(size(chest_ap_outer_mm)); -chest_ap_outer_mm; zeros(size(chest_ap_outer_mm))];
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(thoraxCenter + [0; bodyShift; 0], [0, 0, 0]));
            thorax = MultipleMaterialPhantom([beatingHeart, breathingLung, fatComposite, tissueComposite], ...
                thoraxPose);

            breast_gap_mm = 60;
            breast_radius_mm = 60;
            breast_depth_mm = 125;

            right_breast_center = [breast_radius_mm + 0.5 * breast_gap_mm, 0, 0];
            left_breast_center = [-right_breast_center(1), right_breast_center(2:3)];

            breastParams = struct('radius_mm', breast_radius_mm, 'length_mm', breast_depth_mm);
            breastParamsRight = BreastPhantom.addPoseToParameters(breastParams, right_breast_center, [0, 90, 90]);
            breastParamsLeft = BreastPhantom.addPoseToParameters(breastParams, left_breast_center, [0, 90, 90]);
            breast_right = AnalyticalCylinder3D([], breastParamsRight);
            breast_left = AnalyticalCylinder3D([], breastParamsLeft);

            breastCenter = [0, 0.5 * breast_depth_mm, 0];
            breastParamsBoth = BreastPhantom.addPoseToParameters(breastParams, breastCenter, [0, 0, 0]);

            vesselDiameter_mm = 5;
            vesselRadius_mm = 0.5 * vesselDiameter_mm;
            totalVesselLength_mm = 100;
            breastRollPitchYaw = [0, 90, 90];
            enhancedIntensity = 2.5;
            unenhancedIntensity = 0.3;
            breastVesselVelocity_cm_s = 1/10.66;
            startInjectionTime_s = 30*10.66;
            
            min_max_t = [min(t_s) max(t_s)]

            enhancedParamsHandle = @() localVesselParameters('enhanced', t_s);
            unenhancedParamsHandle = @() localVesselParameters('unenhanced', t_s);

            enhancedCylinder = AnalyticalCylinder3D(enhancedIntensity, enhancedParamsHandle);
            unenhancedCylinder = AnalyticalCylinder3D(unenhancedIntensity, unenhancedParamsHandle);

            % enhancingVessel = MultipleMaterialPhantom([unenhancedCylinder, enhancedCylinder]);
            % enhancingVessel = MultipleMaterialPhantom([unenhancedCylinder]);

            breastIntensity = 0.5;
            leftAndRightBreastTissue = CompositeAnalyticalShape3D([breast_right, breast_left], [unenhancedCylinder, enhancedCylinder], ...
                breastIntensity , breastParamsBoth);

            enhancingVessel = MultipleMaterialPhantom([enhancedCylinder, unenhancedCylinder], breastParamsBoth);

            % enhancingVessel  = CompositeAnalyticalShape3D([unenhancedCylinder, enhancedCylinder], [], ...
            %     breastIntensity , breastParamsBoth);

            obj.setShapes([thorax leftAndRightBreastTissue enhancingVessel]);

            function params = localVesselParameters(segment, t_s)
                elapsedTime_s = max(t_s - startInjectionTime_s,0);
                enhancedLength_mm = min(breastVesselVelocity_cm_s*10 .* elapsedTime_s, totalVesselLength_mm);
                unenhancedLength_mm = max(totalVesselLength_mm - enhancedLength_mm, 0);

                switch lower(segment)
                    case 'enhanced'
                        length_mm = enhancedLength_mm;
                        center_mm = repmat(right_breast_center, numel(length_mm), 1)';
                        center_mm(3,:) = center_mm(3,:) - 0.5 .* unenhancedLength_mm;
                    case 'unenhanced'
                        length_mm = unenhancedLength_mm;
                        center_mm = repmat(right_breast_center, numel(length_mm), 1)';
                        center_mm(3,:) = center_mm(3,:) + 0.5 .* enhancedLength_mm;
                    otherwise
                        error('BreastPhantom:InvalidVesselSegment', ...
                            'Segment must be ''enhanced'' or ''unenhanced''.');
                end

                params = struct('radius_mm', vesselRadius_mm, 'length_mm', length_mm);
                params.pose = BreastPhantom.createPoseStruct(center_mm, breastRollPitchYaw);
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

end
