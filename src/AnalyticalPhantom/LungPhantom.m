classdef LungPhantom < MultipleMaterialPhantom
    % LungPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalLungPhantom.m setup.

    methods
        function obj = LungPhantom()
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

            % A/P blood vessel
            vessel_radius_mm = 2.5;
            total_vessel_length_mm = 100;
            enhancedLength = 10;
            unenhancedLength = total_vessel_length_mm - enhancedLength;
            enhancedCenter = [right_breast_center(1), right_breast_center(2) - 0.5 * unenhancedLength, 0];
            unenhancedCenter = [right_breast_center(1), right_breast_center(2) + 0.5 * enhancedLength, 0];

            vessel_ap_unenhanced = AnalyticalCylinder3D(vessel_radius_mm, unenhancedLength, 0.4, ...
                unenhancedCenter, [0, 90, 90]);
            vessel_ap_enhanced = AnalyticalCylinder3D(vessel_radius_mm, enhancedLength, 2.5, ...
                enhancedCenter, [0, 90, 90]);

            breastRightTissue = CompositeAnalyticalShape3D(breast_right, ...
                [vessel_ap_unenhanced, vessel_ap_enhanced], 0.5, [], []);

            shapes = [fatComposite, tissueComposite, heart, rightLung, leftLung, ...
                breast_left, breastRightTissue, vessel_ap_unenhanced, vessel_ap_enhanced];

            obj@MultipleMaterialPhantom(shapes);
        end
    end
end
