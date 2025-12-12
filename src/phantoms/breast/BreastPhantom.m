classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup. Heart geometry leverages
    %   cardiac_ellipsoid_waveform, which evaluates long time vectors in
    %   ~0.001 Gb, 15,625-sample chunks to keep memory usage minimal.

    properties (Access = private)
        time_s (:,1) double {mustBeFinite}
    end

    methods
        function obj = BreastPhantom(t_s)
            arguments
                t_s (:,1) double {mustBeFinite}
            end

            obj@MultipleMaterialPhantom();

            obj.time_s = t_s(:).';

            %% Create heart
            cardiacOpts = struct('HR_bpm', 70, ...
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
            pulmonaryOpts = struct('f_bpm', 12, ...
                'VT_L', 0.4,...
                'Vres_L', 0.8, ...
                'Vbase_L', 1.5, ...
                'bellyFrac', 0.6, ...
                'inspFrac', 1/3, ...
                'GCS_peak', 0.4);
            lungIntensity = 0.1;
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
            
            thoraxCenter = [zeros(size(chest_ap_outer_mm(:))), -chest_ap_outer_mm(:), zeros(size(chest_ap_outer_mm(:)))]';
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(thoraxCenter + [0; bodyShift; 0], [0, 0, 0]));
            thorax = MultipleMaterialPhantom([beatingHeart, breathingLung, fatComposite, tissueComposite], ...
                thoraxPose);

           
            % 
            % [breastLeft, breastRight, breastCenter, rightBreastCenter] = obj.createBreastGeometry();
            % enhancingVessel = obj.createEnhancingVessel(obj.time_s, rightBreastCenter);
            % 
            % leftAndRightBreastTissue = CompositeAnalyticalShape3D([breastRight, breastLeft], enhancingVessel, ...
            %     0.5, [0, 0, 0], [0, 0, 0]);
            % 
            % bothBreasts = MultipleMaterialPhantom([leftAndRightBreastTissue, enhancingVessel], ...
            %     breastCenter, [0, 0, 0]);

            obj.setShapes([thorax]);
        end
    end

    methods (Access = private)

        function [breastLeft, breastRight, breastCenter, right_breast_center] = createBreastGeometry(~)
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
        end

        function enhancingVessel = createEnhancingVessel(obj, rightBreastCenter)
            vesselRadius_mm = 2.5;
            total_vessel_length_mm = 100;
            rollPitchYaw = [0, 90, 90];

            totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
            V_contrast_mm3 = totalVolume_mm3 * ones(numel(obj.time_s), 1);

            enhancingVessel = EnhancingVessel(obj.time_s.', total_vessel_length_mm, 2.5, 0.4, ...
                vesselRadius_mm, V_contrast_mm3, rightBreastCenter, rollPitchYaw);
        end
    end

    methods (Access = private, Static)
        function params = addPoseToParameters(params, centerVec, rollPitchYaw)
            params.pose = BreastPhantom.createPoseStruct(centerVec, rollPitchYaw);
        end

        function pose = createPoseStruct(centerVec, rollPitchYaw)
            pose = struct('center', struct('x_mm', centerVec(:,1), ...
                'y_mm', centerVec(:,2), ...
                'z_mm', centerVec(:,3)), ...
                'roll_deg', rollPitchYaw(1), ...
                'pitch_deg', rollPitchYaw(2), ...
                'yaw_deg', rollPitchYaw(3));
        end
    end

end
