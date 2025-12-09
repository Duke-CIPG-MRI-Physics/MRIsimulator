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
            % Defaults
            tempCenter = [0, 0, 0]; % will be corrected in future 
            noRotation = [0, 0, 0];

            obj@MultipleMaterialPhantom();

            t_row = t_s(:).';
            obj.time_s = t_row;



            % Heart parameters
            heartOpts = struct('systFrac', 0.35, ...
                'q_ED', 50/27, ...
                'GLS_peak', -0.20, ...
                'GCS_peak', -0.25);
            HR_bpm = 70 * ones(1, numel(t_row));
            EDV_ml = 150 * ones(1, numel(t_row));
            ESV_ml = 75  * ones(1, numel(t_row));

            % Heart
            heart = BeatingHeart(t_row, HR_bpm, EDV_ml, ESV_ml, 1, tempCenter, noRotation, heartOpts);
            heartAFcn = BreastPhantom.buildParameterFunction(t_row, heart.getA());
            heartBFcn = BreastPhantom.buildParameterFunction(t_row, heart.getB());
            heartCFcn = BreastPhantom.buildParameterFunction(t_row, heart.getC());
            BreastPhantom.assignParameterFunctions(heart, struct('a', heartAFcn, ...
                'b', heartBFcn, 'c', heartCFcn));

            % Lung parameters
            f_bpm = 12 * ones(1, numel(t_row));
            VT_L = 0.4 * ones(1, numel(t_row));
            Vres_L = 0.8 * ones(1, numel(t_row));
            Vbase_L = 1.5 * ones(1, numel(t_row));
            bellyFrac = 0.6*ones(1, numel(t_row));
            inspFrac = 0.4 * ones(1, numel(t_row));

            % Lung
            heart_lr_mm = heart.getA(); % Lungs get pushed L/R with cardiac cycle
            heartThickness_mm = 8;
            spacingBetweenLungs = heart_lr_mm + heartThickness_mm;
            breathingLung = BreathingLung(t_row, f_bpm, VT_L, Vres_L, ...
                Vbase_L, bellyFrac, inspFrac, spacingBetweenLungs, ...
                0.1, tempCenter, noRotation);
            lungRadiusFcn = BreastPhantom.buildParameterFunction(t_row, breathingLung.getLungRadiusMm());
            lungHeightFcn = BreastPhantom.buildParameterFunction(t_row, breathingLung.getLungHeightMm());
            spacingBetweenLungsFcn = @(time) heartAFcn(time) + heartThickness_mm;
            BreastPhantom.assignParameterFunctions(breathingLung, struct(...
                'radius', lungRadiusFcn, 'height', lungHeightFcn, ...
                'separation', spacingBetweenLungsFcn));

            figure();
            plot(t_row,heart.getA(),'-r');
            hold on
            plot(t_row,heart.getB(),'--k');
            plot(t_row,heart.getC(),':r');

            plot(t_row,breathingLung.getLungRadiusMm,'-b');
            plot(t_row,breathingLung.getLungHeightMm,'--b');
            legend('Heart A','Heart B','Heart C','Lung A','Lung B')



            % TODO - make this depend on breathing and cardiac motion
            bodyShift = -80;
            tissueGap_lr_mm = 30;
            chestApInnerFcn = @(time) max(heartBFcn(time), lungRadiusFcn(time)) + tissueGap_lr_mm;
            chestLrInnerFcn = @(time) (2*lungRadiusFcn(time) + spacingBetweenLungsFcn(time) + tissueGap_lr_mm);
            phantomDepth_mm = 300;
            bodyCenter = [0 bodyShift 0];
            chest_ap_inner_mm = chestApInnerFcn(t_row);
            chest_lr_inner_mm = chestLrInnerFcn(t_row);
            innerDepthFcn = @(time) 0.9 * phantomDepth_mm + zeros(size(time));
            fat_inner = AnalyticalEllipticalCylinder3D(chest_lr_inner_mm, ...
                chest_ap_inner_mm, innerDepthFcn(t_row), [], tempCenter, noRotation);
            BreastPhantom.assignParameterFunctions(fat_inner, struct(...
                'a', chestLrInnerFcn, 'b', chestApInnerFcn, 'c', innerDepthFcn));

            % Peripheral fat (outer - inner shell)
            fatThickness_mm = 10;
            chestApOuterFcn = @(time) chestApInnerFcn(time) + fatThickness_mm;
            chestLrOuterFcn = @(time) chestLrInnerFcn(time) + fatThickness_mm;
            chest_ap_outer_mm = chestApOuterFcn(t_row);
            chest_lr_outer_mm = chestLrOuterFcn(t_row);
            fat_outer = AnalyticalEllipticalCylinder3D(chest_lr_outer_mm, ...
                chest_ap_outer_mm, innerDepthFcn(t_row), [], tempCenter, noRotation);
            BreastPhantom.assignParameterFunctions(fat_outer, struct(...
                'a', chestLrOuterFcn, 'b', chestApOuterFcn, 'c', innerDepthFcn));

            fatComposite = CompositeAnalyticalShape3D(fat_outer, fat_inner, 2, [], []);
            tissueComposite = CompositeAnalyticalShape3D(fat_inner, [heart, breathingLung], ...
                0.5, tempCenter, noRotation);

            thoraxCenter = [zeros(size(chest_ap_outer_mm(:))), -chest_ap_outer_mm(:), zeros(size(chest_ap_outer_mm(:)))];
            thorax = MultipleMaterialPhantom([heart, breathingLung, fatComposite, tissueComposite],...
                thoraxCenter, noRotation);

            % Breasts
            breast_gap_mm = 60;
            breast_radius_mm = 60;
            breast_depth_mm = 125;
            breastRadiusFcn = @(time) breast_radius_mm + zeros(size(time));
            breastDepthFcn = @(time) breast_depth_mm + zeros(size(time));
            right_breast_center = [breast_radius_mm + 0.5 * breast_gap_mm, ...
                0, 0];

            breast_right = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                right_breast_center, [0, 90, 90]);
            BreastPhantom.assignParameterFunctions(breast_right, struct(...
                'radius', breastRadiusFcn, 'height', breastDepthFcn));

            left_breast_center = [-right_breast_center(1), right_breast_center(2:3)];
            breast_left = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                left_breast_center, [0, 90, 90]);
            BreastPhantom.assignParameterFunctions(breast_left, struct(...
                'radius', breastRadiusFcn, 'height', breastDepthFcn));

            % A/P blood vessel with contrast wash-in
            vesselRadius_mm = 2.5;
            total_vessel_length_mm = 100;
            rollPitchYaw = [0, 90, 90];

            % Calculate contrast wash-in
            ContrastStartTime = 0.25 * t_s(end);
            ContrastEndTime   = 0.75 * t_s(end);
            totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
            
            % Option 1: Contrast already present
            V_contrast_mm3 = totalVolume_mm3*ones(numel(t_s), 1);
            
            % Option 2: Contrast linear washin
            % V_contrast_mm3 = zeros(numel(t_s), 1);
            % midRamp = t_s >= ContrastStartTime & t_s <= ContrastEndTime;
            % V_contrast_mm3(midRamp) = totalVolume_mm3 * (t_s(midRamp) - ContrastStartTime) ./ (ContrastEndTime - ContrastStartTime);
            % V_contrast_mm3(t_s > ContrastEndTime) = totalVolume_mm3;
            
            V_contrast_mm3 = V_contrast_mm3(:);
            if numel(V_contrast_mm3) ~= numel(t_row)
                error('BreastPhantom:ContrastSizeMismatch', ...
                    'V_contrast_mm3 must match the length of t_s.');
            end
            
            vesselRadiusFcn = @(time) vesselRadius_mm + zeros(size(time));
            vesselLengthFcn = @(time) total_vessel_length_mm + zeros(size(time));
            enhancingVessel = EnhancingVessel(t_row.', total_vessel_length_mm, 2.5, 0.4, ...
                vesselRadius_mm, V_contrast_mm3, right_breast_center, rollPitchYaw);
            BreastPhantom.assignParameterFunctions(enhancingVessel, struct(...
                'radius', vesselRadiusFcn, 'length', vesselLengthFcn));

            leftAndRightBreastTissue = CompositeAnalyticalShape3D([breast_right, breast_left], enhancingVessel, ...
                0.5, tempCenter, noRotation);

            breastCenter = [0, 0.5*breast_depth_mm, 0];
            bothBreasts = MultipleMaterialPhantom([leftAndRightBreastTissue, enhancingVessel],...
                breastCenter, noRotation);
            
            obj.setShapes([thorax, bothBreasts]);
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
