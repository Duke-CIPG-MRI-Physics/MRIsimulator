classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup. Heart geometry leverages
    %   cardiac_ellipsoid_waveform, which evaluates long time vectors in
    %   ~0.001 Gb, 15,625-sample chunks to keep memory usage minimal.

    properties (Access = protected)
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

            % Convenience arrays
            centered = [0, 0, 0];
            notRotated = [0, 0, 0];

            %cardiacParams = @()BreastPhantom.addPoseToParameters( ...
            %     cardiac_ellipsoid_waveform(t_s, cardiacOpts), centered, notRotated);
            cardiacParamsInit = BreastPhantom.addPoseToParameters( ...
                cardiac_ellipsoid_waveform(0, params.cardiacOpts), centered, notRotated);
            obj.beatingHeart = AnalyticalEllipsoid3D(params.heartIntensity, cardiacParamsInit);

            %% Create lungs
            % TODO make this a function handle so its not directly
            % evaluated yet.
            lungShapeParameters = BreastPhantom.addPoseToParameters(struct(), centered, notRotated);

            [lung_radius_mm, lung_heigh_mm] = lung_ellipsoid_waveform(0, params.pulmonaryOpts);
            % [R_lung_mm, H_Lung_mm] = lung_ellipsoid_waveform(t_s, pulmonaryOpts);
            
            lungSeparation_mm = obj.beatingHeart.getShapeParameters.a_mm + params.heartWallThickness_mm; 
            lungPosition_mm = lung_radius_mm + lungSeparation_mm;

            lungShapeParams = struct('a_mm', lung_radius_mm, 'b_mm', lung_radius_mm, 'c_mm', lung_heigh_mm)
            rightLungParams = lungShapeParams;
            rightLungCenter= struct('x_mm', lungPosition_mm, ...
                'y_mm', 0, ...
                'z_mm', 0);
            rightLungParams.pose = struct('center', rightLungCenter,...
                'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);

            leftLungParams = lungShapeParams;
            leftLungCenter= struct('x_mm', -lungPosition_mm, ...
                'y_mm', 0, ...
                'z_mm', 0);
            leftLungParams.pose = struct('center', leftLungCenter,...
                'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);

            obj.rightLung = AnalyticalEllipsoid3D(params.lungIntensity, rightLungParams);
            obj.leftLung = AnalyticalEllipsoid3D(params.lungIntensity, leftLungParams);
            

            %% Create thorax
            heart_ap_mm = obj.beatingHeart.getShapeParameters.b_mm;
            
            chest_ap_inner_mm = max(heart_ap_mm, lung_radius_mm) + params.tissueGap_lr_mm;
            chest_lr_inner_mm = 2 * lung_radius_mm + lungSeparation_mm + params.tissueGap_lr_mm;
            
            fatInnerParams = struct('a_mm', chest_lr_inner_mm, ...
                'b_mm', chest_ap_inner_mm, 'length_mm', params.phantomDepth_mm);
            fatInnerParams = BreastPhantom.addPoseToParameters(fatInnerParams, [0, 0, 0], [0, 0, 0]);
            obj.fatInner = AnalyticalEllipticalCylinder3D([], fatInnerParams);

            fatThickness_mm = 10;
            chest_ap_outer_mm = chest_ap_inner_mm + fatThickness_mm;
            chest_lr_outer_mm = chest_lr_inner_mm + fatThickness_mm;
            fatOuterParams = struct('a_mm', chest_lr_outer_mm, ...
                'b_mm', chest_ap_outer_mm, 'length_mm', params.phantomDepth_mm);
            fatOuterParams = BreastPhantom.addPoseToParameters(fatOuterParams, [0, 0, 0], [0, 0, 0]);
            obj.fatOuter = AnalyticalEllipticalCylinder3D([], fatOuterParams);

            obj.fatComposite = CompositeAnalyticalShape3D(obj.fatOuter, obj.fatInner, params.fatIntensity);
            obj.tissueComposite = CompositeAnalyticalShape3D(obj.fatInner, [obj.beatingHeart, obj.leftLung, obj.rightLung], ...
               params.tissueIntensity);
            
            thoraxCenter = [zeros(size(chest_ap_outer_mm)); -chest_ap_outer_mm; zeros(size(chest_ap_outer_mm))];
            thoraxPose = struct('pose', BreastPhantom.createPoseStruct(thoraxCenter, [0, 0, 0]));
            obj.thorax = MultipleMaterialPhantom([obj.beatingHeart, obj.leftLung, obj.rightLung, obj.fatComposite, obj.tissueComposite], ...
                thoraxPose);

           
            right_breast_center = [params.breast_radius_mm + 0.5 * params.breast_gap_mm, 0, 0];
            left_breast_center = [-right_breast_center(1), right_breast_center(2:3)];

            breastParams = struct('radius_mm', params.breast_radius_mm, 'length_mm', params.breast_depth_mm);
            breastParamsRight = BreastPhantom.addPoseToParameters(breastParams, right_breast_center, [0, 90, 90]);
            breastParamsLeft = BreastPhantom.addPoseToParameters(breastParams, left_breast_center, [0, 90, 90]);
            obj.breastRight = AnalyticalCylinder3D([], breastParamsRight);
            obj.breastLeft = AnalyticalCylinder3D([], breastParamsLeft);

            breastCenter = [0, 0.5 * params.breast_depth_mm, 0];
            breastParamsBoth = BreastPhantom.addPoseToParameters(breastParams, breastCenter, [0, 0, 0]);
            
            % min_max_t = [min(t_s) max(t_s)]
            min_max_t = [0 0];

            % enhancedParamsHandle = @() localVesselParameters('enhanced', t_s);
            % unenhancedParamsHandle = @() localVesselParameters('unenhanced', t_s);
            enhancedParamsHandle = @() localVesselParameters('enhanced', 0, params);
            unenhancedParamsHandle = @() localVesselParameters('unenhanced', 0, params);

            obj.enhancedCylinder = AnalyticalCylinder3D(params.enhancedIntensity, enhancedParamsHandle);
            obj.unenhancedCylinder = AnalyticalCylinder3D(params.unenhancedIntensity, unenhancedParamsHandle);

            % enhancingVessel = MultipleMaterialPhantom([unenhancedCylinder, enhancedCylinder]);
            % enhancingVessel = MultipleMaterialPhantom([unenhancedCylinder]);

            
            obj.leftAndRightBreastTissue = CompositeAnalyticalShape3D([obj.breastRight, obj.breastLeft], [obj.unenhancedCylinder, obj.enhancedCylinder], ...
                params.breastIntensity , breastParamsBoth);

            obj.enhancingVessel = MultipleMaterialPhantom([obj.enhancedCylinder, obj.unenhancedCylinder], breastParamsBoth);

            % enhancingVessel  = CompositeAnalyticalShape3D([unenhancedCylinder, enhancedCylinder], [], ...
            %     breastIntensity , breastParamsBoth);

            obj.setShapes([obj.thorax obj.leftAndRightBreastTissue obj.enhancingVessel]);

            function vesselParams = localVesselParameters(segment, t_s, params)
                elapsedTime_s = max(t_s - params.startInjectionTime_s,0);
                enhancedLength_mm = min(params.breastVesselVelocity_cm_s*10 .* elapsedTime_s, params.totalVesselLength_mm);
                unenhancedLength_mm = max(params.totalVesselLength_mm - enhancedLength_mm, 0);

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

                vesselParams = struct('radius_mm', params.vesselRadius_mm, 'length_mm', length_mm);
                vesselParams.pose = BreastPhantom.createPoseStruct(center_mm, params.breastRollPitchYaw);
            end
        end

        function S = kspaceAtTime(obj, kx, ky, kz, t_s)
            S = obj.kspace(kx, ky, kz);
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
