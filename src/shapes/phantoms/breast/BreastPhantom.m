classdef BreastPhantom < MultiIntensityShapeGroup3D
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup.

    methods
        function obj = BreastPhantom(provider)
            arguments
                provider (1,1) PhantomContext
            end
            % Defaults
            tempCenter = [0, 0, 0]; % will be corrected in future 
            noRotation = [0, 0, 0];

            obj@MultiIntensityShapeGroup3D();

            t_row = provider.getTime();

            % Heart
            heart = BeatingHeart(provider, 1, tempCenter, noRotation);

            % Lung
            [heartA, heartB, ~] = provider.getHeartRadiiMm(); % Lungs get pushed L/R with cardiac cycle
            heartThickness_mm = 8;
            spacingBetweenLungs = @(t) heartA(t) + heartThickness_mm;
            breathingLung = BreathingLung(provider, spacingBetweenLungs, ...
                0.1, tempCenter, noRotation);
            


            % TODO - make this depend on breathing and cardiac motion
            bodyShift = -80;
            heart_ap_mm = heartB(t_row);
            lungRadius = breathingLung.getLungRadiusMm(t_row);
            tissueGap_lr_mm = 20;
            tissueGap_lr_mm = 30;
            chest_ap_inner_mm = (max(heart_ap_mm,lungRadius) + tissueGap_lr_mm);
            spacingBetweenLungs_mm = spacingBetweenLungs(t_row);
            chest_lr_inner_mm = (2*lungRadius + spacingBetweenLungs_mm + tissueGap_lr_mm);
            phantomDepth_mm = 300;
            bodyCenter = [0 bodyShift 0];
            fat_inner = AnalyticalEllipticalCylinder3D(chest_lr_inner_mm, ...
                chest_ap_inner_mm, 0.9 * phantomDepth_mm, [], tempCenter, noRotation);

            % Peripheral fat (outer - inner shell)
            fatThickness_mm = 10;
            chest_ap_outer_mm = chest_ap_inner_mm + fatThickness_mm;
            chest_lr_outer_mm = chest_lr_inner_mm + fatThickness_mm;
            fat_outer = AnalyticalEllipticalCylinder3D(chest_lr_outer_mm, ...
                chest_ap_outer_mm, 0.9 * phantomDepth_mm, [], tempCenter, noRotation);

            fatComposite = SharedIntensityShapeGroup3D(fat_outer, fat_inner, 2, [], []);
            tissueComposite = SharedIntensityShapeGroup3D(fat_inner, [heart, breathingLung], ...
                0.5, tempCenter, noRotation);

            thoraxCenter = [zeros(size(chest_ap_outer_mm(:))), -chest_ap_outer_mm(:), zeros(size(chest_ap_outer_mm(:)))];
            thorax = MultiIntensityShapeGroup3D([heart, breathingLung, fatComposite, tissueComposite],...
                thoraxCenter, noRotation);

            % Breasts
            breast_gap_mm = 60;
            breast_radius_mm = 60;
            breast_depth_mm = 125;
            right_breast_center = [breast_radius_mm + 0.5 * breast_gap_mm, ...
                0, 0];

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
            ContrastStartTime = 0.25 * t_row(end);
            ContrastEndTime   = 0.75 * t_row(end);
            totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;

            % Option 1: Contrast already present
            V_contrast_mm3 = provider.getVesselContrastMm3();
            V_contrast_mm3 = V_contrast_mm3(t_row).';
            
            % Option 2: Contrast linear washin
            % V_contrast_mm3 = zeros(numel(t_row), 1);
            % midRamp = t_row >= ContrastStartTime & t_row <= ContrastEndTime;
            % V_contrast_mm3(midRamp) = totalVolume_mm3 * (t_row(midRamp) - ContrastStartTime) ./ (ContrastEndTime - ContrastStartTime);
            % V_contrast_mm3(t_row > ContrastEndTime) = totalVolume_mm3;

            V_contrast_mm3 = V_contrast_mm3(:);
            if numel(V_contrast_mm3) ~= numel(t_row)
                error('BreastPhantom:ContrastSizeMismatch', ...
                    'V_contrast_mm3 must match the length of the time vector.');
            end
            
            enhancingVessel = EnhancingVessel(t_row.', total_vessel_length_mm, 2.5, 0.4, ...
                vesselRadius_mm, V_contrast_mm3, right_breast_center, rollPitchYaw);

            leftAndRightBreastTissue = SharedIntensityShapeGroup3D([breast_right, breast_left], enhancingVessel, ...
                0.5, tempCenter, noRotation);

            breastCenter = [0, 0.5*breast_depth_mm, 0];
            bothBreasts = MultiIntensityShapeGroup3D([leftAndRightBreastTissue, enhancingVessel],...
                breastCenter, noRotation);
            
            obj.setShapes([thorax, bothBreasts]);
        end
    end
end
