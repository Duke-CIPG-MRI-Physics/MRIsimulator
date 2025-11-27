classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup.

    methods
        function obj = BreastPhantom(t_s, V_contrast_mm3, vesselRadius_mm)
            arguments
                t_s (:,1) double {mustBeFinite}
                V_contrast_mm3 (:,1) double {mustBeFinite, mustBeNonnegative} = [];
                vesselRadius_mm double {mustBePositive} = 2.5;
            end

            bodyShift = -80;

            % Heart
            heart_a_mm = 50;
            heart_b_mm = 27;
            heart_c_mm = 30;
            heart_center = [0 bodyShift 0];
            heart = AnalyticalEllipsoid3D(heart_a_mm, heart_b_mm, heart_c_mm, 1, ...
                heart_center, [0, -65, 70]);

            % Lungs
            lung_a_mm = 140;
            lung_b_mm = 50;
            lung_c_mm = 50;
            lungSeparation = max(lung_b_mm, lung_c_mm) + max(heart_b_mm, heart_c_mm) + 2;

            center_R_mm = [lungSeparation, bodyShift, 0];
            rightLung = AnalyticalEllipsoid3D(lung_a_mm, lung_b_mm, lung_c_mm, 0.1, ...
                center_R_mm, [0, 95, 0]);

            center_L_mm = [-lungSeparation, bodyShift, 0];
            leftLung = AnalyticalEllipsoid3D(lung_a_mm, lung_b_mm, lung_c_mm, 0.1, ...
                center_L_mm, [0, 85, 0]);

            % Peripheral fat (outer - inner shell)
            fatThickness_mm = 10;
            tissueThickness_mm = 10;
            patientThickness_outer_mm = 1.85 * (max(lung_b_mm, lung_c_mm) + fatThickness_mm);
            patientWidth_outer_mm = 1.5 * (lungSeparation + max(lung_b_mm, lung_c_mm) + ...
                fatThickness_mm + tissueThickness_mm) + 2;

            bodyCenter = [0 bodyShift 0];
            phantomDepth_mm = 400;
            fat_outer = AnalyticalEllipticalCylinder3D(patientWidth_outer_mm, ...
                patientThickness_outer_mm, 0.9 * phantomDepth_mm, [], bodyCenter, [0, 0, 0]);

            patientThickness_inner_mm = patientThickness_outer_mm - 2 * fatThickness_mm;
            patientWidth_inner_mm = patientWidth_outer_mm - 2 * fatThickness_mm;
            fat_inner = AnalyticalEllipticalCylinder3D(patientWidth_inner_mm, ...
                patientThickness_inner_mm, 0.9 * phantomDepth_mm, [], bodyCenter, [0, 0, 0]);

            fatComposite = CompositeAnalyticalShape3D(fat_outer, fat_inner, 2, [], []);
            tissueComposite = CompositeAnalyticalShape3D(fat_inner, [heart, rightLung, leftLung], ...
                0.5, [], []);

            % Breasts
            breast_gap_mm = 50;
            breast_radius_mm = 65;
            breast_depth_mm = 200;
            right_breast_center = [breast_radius_mm + 0.5 * breast_gap_mm, ...
                bodyShift + 0.5 * breast_depth_mm + patientThickness_outer_mm, 0];

            breast_right = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, [], ...
                right_breast_center, [0, 90, 90]);
            breast_left = AnalyticalCylinder3D(breast_radius_mm, breast_depth_mm, 0.5, ...
                [-right_breast_center(1), right_breast_center(2:3)], [0, 90, 90]);

            % A/P blood vessel with contrast wash-in
            total_vessel_length_mm = 100;
            rollPitchYaw = [0, 90, 90];

            if isempty(V_contrast_mm3)
                totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
                V_contrast_mm3 = linspace(0, totalVolume_mm3, numel(t_s)).';
            else
                V_contrast_mm3 = V_contrast_mm3(:);
                if numel(V_contrast_mm3) ~= numel(t_s)
                    error('BreastPhantom:ContrastSizeMismatch', ...
                        'V_contrast_mm3 must match the length of t_s.');
                end
            end

            enhancingVessel = EnhancingVessel(t_s, total_vessel_length_mm, 2.5, 0.4, ...
                vesselRadius_mm, V_contrast_mm3, right_breast_center, rollPitchYaw);
            enhancingVessel.updateTimeArray(t_s(end), V_contrast_mm3(end), vesselRadius_mm);
            [enhancedSegment, unenhancedSegment] = enhancingVessel.getVessels();

            % Ensure vessel segments are provided as a row vector of
            % AnalyticalShape3D objects to satisfy CompositeAnalyticalShape3D
            % argument validation even if upstream callers supply column
            % vectors.
            vesselSegments = [unenhancedSegment(:).', enhancedSegment(:).'];

            breastRightTissue = CompositeAnalyticalShape3D(breast_right, ...
                vesselSegments, 0.5, [], []);

            shapes = [fatComposite, tissueComposite, heart, rightLung, leftLung, ...
                breast_left, breastRightTissue, unenhancedSegment, enhancedSegment];

            obj@MultipleMaterialPhantom(shapes);
        end
    end
end
