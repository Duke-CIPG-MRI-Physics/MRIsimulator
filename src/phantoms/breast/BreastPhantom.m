classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup.

    properties (Access = private)
        time_s (:,1) double {mustBeFinite}
    end

    methods
        function obj = BreastPhantom(t_s)
            arguments
                t_s (:,1) double {mustBeFinite}
            end

            obj@MultipleMaterialPhantom();

            t_row = t_s(:).';
            obj.time_s = t_row;

            heart = obj.createHeart(t_row);
            breathingLung = obj.createBreathingLung(t_row, heart);
            thorax = obj.createThorax(heart, breathingLung);

            [breastLeft, breastRight, breastCenter, rightBreastCenter] = obj.createBreastGeometry();
            enhancingVessel = obj.createEnhancingVessel(t_row, rightBreastCenter);

            leftAndRightBreastTissue = CompositeAnalyticalShape3D([breastRight, breastLeft], enhancingVessel, ...
                0.5, [0, 0, 0], [0, 0, 0]);

            bothBreasts = MultipleMaterialPhantom([leftAndRightBreastTissue, enhancingVessel], ...
                breastCenter, [0, 0, 0]);

            obj.setShapes([thorax, bothBreasts]);
        end
    end

    methods (Access = private)
        function heart = createHeart(~, t_row)
            heartOpts = struct('systFrac', 0.35, ...
                'q_ED', 50/27, ...
                'GLS_peak', -0.20, ...
                'GCS_peak', -0.25);
            HR_bpm = 70 * ones(1, numel(t_row));
            EDV_ml = 150 * ones(1, numel(t_row));
            ESV_ml = 75  * ones(1, numel(t_row));

            heart = BeatingHeart(t_row, HR_bpm, EDV_ml, ESV_ml, 1, [0, 0, 0], [0, 0, 0], heartOpts);
        end

        function breathingLung = createBreathingLung(~, t_row, heart)
            f_bpm = 12 * ones(1, numel(t_row));
            VT_L = 0.4 * ones(1, numel(t_row));
            Vres_L = 0.8 * ones(1, numel(t_row));
            Vbase_L = 1.5 * ones(1, numel(t_row));
            bellyFrac = 0.6 * ones(1, numel(t_row));
            inspFrac = 0.4 * ones(1, numel(t_row));

            heart_lr_mm = heart.getA();
            heartThickness_mm = 8;
            spacingBetweenLungs = heart_lr_mm + heartThickness_mm;
            breathingLung = BreathingLung(t_row, f_bpm, VT_L, Vres_L, ...
                Vbase_L, bellyFrac, inspFrac, spacingBetweenLungs, ...
                0.1, [0, 0, 0], [0, 0, 0]);
        end

        function thorax = createThorax(~, heart, breathingLung)
            bodyShift = -80;
            heart_ap_mm = heart.getB();
            lungRadius = breathingLung.getLungRadiusMm();
            tissueGap_lr_mm = 30;
            heart_lr_mm = heart.getA();
            heartThickness_mm = 8;
            spacingBetweenLungs = heart_lr_mm + heartThickness_mm;

            chest_ap_inner_mm = max(heart_ap_mm, lungRadius) + tissueGap_lr_mm;
            chest_lr_inner_mm = 2 * lungRadius + spacingBetweenLungs + tissueGap_lr_mm;
            phantomDepth_mm = 300;

            fat_inner = AnalyticalEllipticalCylinder3D(chest_lr_inner_mm, ...
                chest_ap_inner_mm, 0.9 * phantomDepth_mm, [], [0, 0, 0], [0, 0, 0]);

            fatThickness_mm = 10;
            chestApOuterFcn = @(time) chestApInnerFcn(time) + fatThickness_mm;
            chestLrOuterFcn = @(time) chestLrInnerFcn(time) + fatThickness_mm;
            chest_ap_outer_mm = chestApOuterFcn(t_row);
            chest_lr_outer_mm = chestLrOuterFcn(t_row);
            fat_outer = AnalyticalEllipticalCylinder3D(chest_lr_outer_mm, ...
                chest_ap_outer_mm, 0.9 * phantomDepth_mm, [], [0, 0, 0], [0, 0, 0]);

            fatComposite = CompositeAnalyticalShape3D(fat_outer, fat_inner, 2, [], []);
            tissueComposite = CompositeAnalyticalShape3D(fat_inner, [heart, breathingLung], ...
                0.5, [0, 0, 0], [0, 0, 0]);

            thoraxCenter = [zeros(size(chest_ap_outer_mm(:))), -chest_ap_outer_mm(:), zeros(size(chest_ap_outer_mm(:)))];
            thorax = MultipleMaterialPhantom([heart, breathingLung, fatComposite, tissueComposite], ...
                thoraxCenter + [0, bodyShift, 0], [0, 0, 0]);
        end

        function [breastLeft, breastRight, breastCenter, right_breast_center] = createBreastGeometry(~)
            breast_gap_mm = 60;
            breast_radius_mm = 60;
            breast_depth_mm = 125;

            right_breast_center = [breast_radius_mm + 0.5 * breast_gap_mm, 0, 0];
            left_breast_center = [-right_breast_center(1), right_breast_center(2:3)];

            breast_right = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                right_breast_center, [0, 90, 90]);
            breast_left = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                left_breast_center, [0, 90, 90]);
            BreastPhantom.assignParameterFunctions(breast_left, struct(...
                'radius', breastRadiusFcn, 'height', breastDepthFcn));

            breastCenter = [0, 0.5 * breast_depth_mm, 0];
        end

        function enhancingVessel = createEnhancingVessel(~, t_row, rightBreastCenter)
            vesselRadius_mm = 2.5;
            total_vessel_length_mm = 100;
            rollPitchYaw = [0, 90, 90];

            totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
            V_contrast_mm3 = totalVolume_mm3 * ones(numel(t_row), 1);

            enhancingVessel = EnhancingVessel(t_row.', total_vessel_length_mm, 2.5, 0.4, ...
                vesselRadius_mm, V_contrast_mm3, rightBreastCenter, rollPitchYaw);
        end
    end

    methods (Static, Access = private)
        function parameterFcn = buildParameterFunction(timeSamples, parameterValues)
            arguments
                timeSamples (1,:) double {mustBeFinite}
                parameterValues (1,:) double {mustBeFinite}
            end

            parameterFcn = @(queryTime) interp1(timeSamples(:), parameterValues(:), ...
                queryTime, 'linear', 'extrap');
        end

        function assignParameterFunctions(shape, parameterFunctions)
            if isempty(parameterFunctions)
                return;
            end

            if ismethod(shape, 'setShapeParameterFunctions')
                shape.setShapeParameterFunctions(parameterFunctions);
            elseif isprop(shape, 'shapeParameterFunctions')
                shape.shapeParameterFunctions = parameterFunctions;
            elseif isa(shape, 'dynamicprops')
                newProp = addprop(shape, 'shapeParameterFunctions'); %#ok<NASGU>
                shape.shapeParameterFunctions = parameterFunctions;
            end
        end
    end
end
