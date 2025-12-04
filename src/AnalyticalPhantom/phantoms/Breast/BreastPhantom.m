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

            % Defaults
            tempCenter = [0, 0, 0]; % will be corrected in future 
            noRotation = [0, 0, 0];

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

            % Lung parameters
            f_bpm = 12 * ones(1, numel(t_row));
            VT_L = 0.4 * ones(1, numel(t_row));
            Vres_L = 0.8 * ones(1, numel(t_row));
            Vbase_L = 1.5 * ones(1, numel(t_row));
            bellyFrac = zeros(1, numel(t_row));
            inspFrac = 0.4 * ones(1, numel(t_row));

            % Lung
            heart_lr_mm = heart.getA(); % Lungs get pushed L/R with cardiac cycle
            heartThickness_mm = 8;
            spacingBetweenLungs = heart_lr_mm + heartThickness_mm;
            breathingLung = BreathingLung(t_row, f_bpm, VT_L, Vres_L, ...
                Vbase_L, bellyFrac, inspFrac, spacingBetweenLungs, ...
                0.1, tempCenter, noRotation);

            % TODO - make this depend on breathing and cardiac motion
            bodyShift = -80;
            heart_ap_mm = heart.getB();

            

            lungRadius = breathingLung.getLungRadiusMm();
            maxLungSize = breathingLung.getMaxLungSizeMm();
            lungSeparation = breathingLung.getLungSeparationMm();

            % Peripheral fat (outer - inner shell)
            fatThickness_mm = 10;
            tissueThickness_mm = 10;
            patientThickness_outer_mm = 1.85 * (lungRadius + fatThickness_mm);
            patientWidth_outer_mm = lungSeparation + lungRadius + ...
                fatThickness_mm + tissueThickness_mm + 1;

            bodyCenter = [0 bodyShift 0];
            phantomDepth_mm = 400;
            fat_outer = AnalyticalEllipticalCylinder3D(patientWidth_outer_mm, ...
                patientThickness_outer_mm, 0.9 * phantomDepth_mm, [], tempCenter, noRotation);

            patientThickness_inner_mm = patientThickness_outer_mm - 2 * fatThickness_mm;
            patientWidth_inner_mm = patientWidth_outer_mm - 2 * fatThickness_mm;
            fat_inner = AnalyticalEllipticalCylinder3D(patientWidth_inner_mm, ...
                patientThickness_inner_mm, 0.9 * phantomDepth_mm, [], tempCenter, noRotation);

            fatComposite = CompositeAnalyticalShape3D(fat_outer, fat_inner, 2, [], []);
            tissueComposite = CompositeAnalyticalShape3D(fat_inner, [heart, breathingLung], ...
                0.5, tempCenter, noRotation);

            % Breasts
            breast_gap_mm = 50;
            breast_radius_mm = 65;
            breast_depth_mm = 200;
            right_breast_center = [breast_radius_mm + 0.5 * breast_gap_mm, ...
                bodyShift + 0.5 * breast_depth_mm + patientThickness_outer_mm, 0];

            breast_right = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                right_breast_center, [0, 90, 90]);

            left_breast_center = [-right_breast_center(1), right_breast_center(2:3)];
            breast_left = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                left_breast_center, [0, 90, 90]);

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
            
            enhancingVessel = EnhancingVessel(t_row.', total_vessel_length_mm, 2.5, 0.4, ...
                vesselRadius_mm, V_contrast_mm3, right_breast_center, rollPitchYaw);

            rightBreastComposite = CompositeAnalyticalShape3D(breast_right, enhancingVessel, ...
                0.5, tempCenter, noRotation);

            % obj.setShapes([breathingLung]);
            obj.setShapes([heart, breathingLung, fatComposite, tissueComposite, breast_left, rightBreastComposite, enhancingVessel]);
        end
    end
end
