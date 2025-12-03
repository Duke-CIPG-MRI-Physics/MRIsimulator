classdef BreastPhantom < MultipleMaterialPhantom
    % BreastPhantom
    %   Preconfigured collection of AnalyticalShape3D objects approximating
    %   a thoracic slice with lungs, heart, peripheral fat, breasts, and a
    %   simple vessel. Geometry and intensities follow the original
    %   demo_analyticalBreastPhantom.m setup.

    properties (Access = private)
        time_s (:,1) double {mustBeFinite}
        contrastVolume_mm3 (:,1) double {mustBeFinite, mustBeNonnegative}
        vesselRadius_mm (:,1) double {mustBePositive}
        enhancingVessel (1,1) EnhancingVessel = EnhancingVessel.empty
        totalVesselLength_mm double {mustBePositive} = 100;
    end

    methods
        function obj = BreastPhantom(t_s, V_contrast_mm3, vesselRadius_mm)
            arguments
                t_s (:,1) double {mustBeFinite}
                V_contrast_mm3 (:,1) double {mustBeFinite, mustBeNonnegative} = [];
                vesselRadius_mm double {mustBePositive} = 2.5;
            end

            obj@MultipleMaterialPhantom();

            bodyShift = -80;

            t_row = t_s(:).';

            heartOpts = struct('systFrac', 0.35, ...
                'q_ED', 50/27, ...
                'GLS_peak', -0.20, ...
                'GCS_peak', -0.25);

            HR_bpm = 70 * ones(1, numel(t_row));
            EDV_ml = 150 * ones(1, numel(t_row));
            ESV_ml = 75  * ones(1, numel(t_row));

            % Heart
            heart_center = [0 bodyShift 0];
            heart = BeatingHeart(t_row, HR_bpm, EDV_ml, ESV_ml, 1, heart_center, ...
                [0, -65, 70], heartOpts);

            heart_b_mm = heart.getB();
            heart_c_mm = heart.getC();
            heart_b_max_mm = max(heart_b_mm(:));
            heart_c_max_mm = max(heart_c_mm(:));

            maxHeartDim_mm = max(heart_b_max_mm, heart_c_max_mm);

            obj.time_s = t_s(:);
            obj.vesselRadius_mm = obj.ensureRadiusVector(vesselRadius_mm, numel(obj.time_s));

            % Lungs
            f_bpm = 12 * ones(1, numel(t_row));
            VT_L = 0.4 * ones(1, numel(t_row));
            Vres_L = 0.8 * ones(1, numel(t_row));
            Vbase_L = 1.5 * ones(1, numel(t_row));
            bellyFrac = zeros(1, numel(t_row));
            inspFrac = 0.4 * ones(1, numel(t_row));

            breathingLung = BreathingLung(t_row, f_bpm, VT_L, Vres_L, ...
                Vbase_L, bellyFrac, inspFrac, maxHeartDim_mm, 0.1, [0, bodyShift, 0]);

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
                patientThickness_outer_mm, 0.9 * phantomDepth_mm, [], bodyCenter, [0, 0, 0]);

            patientThickness_inner_mm = patientThickness_outer_mm - 2 * fatThickness_mm;
            patientWidth_inner_mm = patientWidth_outer_mm - 2 * fatThickness_mm;
            fat_inner = AnalyticalEllipticalCylinder3D(patientWidth_inner_mm, ...
                patientThickness_inner_mm, 0.9 * phantomDepth_mm, [], bodyCenter, [0, 0, 0]);

            fatComposite = CompositeAnalyticalShape3D(fat_outer, fat_inner, 2, [], []);
            tissueComposite = CompositeAnalyticalShape3D(fat_inner, [heart, breathingLung], ...
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
            rollPitchYaw = [0, 90, 90];

            obj.contrastVolume_mm3 = obj.ensureContrastCurve(V_contrast_mm3, ...
                numel(obj.time_s));

            obj.enhancingVessel = EnhancingVessel(obj.time_s, obj.totalVesselLength_mm, 2.5, 0.4, ...
                obj.vesselRadius_mm, obj.contrastVolume_mm3, right_breast_center, rollPitchYaw);

            breastRightTissue = CompositeAnalyticalShape3D(breast_right, ...
                obj.enhancingVessel, 0.5, [], []);

            shapes = [fatComposite, tissueComposite, heart, breathingLung, ...
                breast_left, breastRightTissue, obj.enhancingVessel];

            obj@MultipleMaterialPhantom(shapes);
        end

        function t_s = getTimeVector(obj)
            % getTimeVector
            %   Return the time base used when constructing the phantom.
            t_s = obj.time_s;
        end

        function contrastCurve = getContrastVolumeCurve(obj)
            % getContrastVolumeCurve
            %   Return the contrast volume curve driving the enhancing
            %   vessel.
            contrastCurve = obj.contrastVolume_mm3;
        end

        function vesselRadius_mm = getVesselRadius(obj)
            % getVesselRadius
            %   Return the vessel radius used for the enhancing vessel.
            vesselRadius_mm = obj.vesselRadius_mm;
        end

        function enhancingVessel = getEnhancingVessel(obj)
            % getEnhancingVessel
            %   Return the EnhancingVessel instance embedded in the phantom.
            enhancingVessel = obj.enhancingVessel;
        end

        function totalLength_mm = getTotalVesselLength(obj)
            % getTotalVesselLength
            %   Return the total vessel length used in the enhancing vessel model.
            totalLength_mm = obj.totalVesselLength_mm;
        end
    end

    methods (Access = private)
        function radiusVec = ensureRadiusVector(~, radiusInput, targetLength)
            if isscalar(radiusInput)
                radiusVec = radiusInput * ones(targetLength, 1);
            else
                if numel(radiusInput) ~= targetLength
                    error('BreastPhantom:RadiusSizeMismatch', ...
                        'vesselRadius_mm must be scalar or match length of t_s.');
                end
                radiusVec = radiusInput(:);
            end
        end

        function contrastCurve = ensureContrastCurve(obj, contrastInput, targetLength)
            if isempty(contrastInput)
                contrastCurve = obj.defaultContrastCurve();
            else
                contrastCurve = contrastInput(:);
                if numel(contrastCurve) ~= targetLength
                    error('BreastPhantom:ContrastSizeMismatch', ...
                        'V_contrast_mm3 must match the length of t_s.');
                end
            end
        end

        function contrastCurve = defaultContrastCurve(obj)
            startTime = 0.25 * obj.time_s(end);
            endTime = 0.75 * obj.time_s(end);

            contrastCurve = zeros(numel(obj.time_s), 1);

            totalVolume_mm3 = pi .* obj.vesselRadius_mm.^2 .* obj.totalVesselLength_mm;

            midRamp = obj.time_s >= startTime & obj.time_s <= endTime;
            contrastCurve(midRamp) = totalVolume_mm3 .* ...
                (obj.time_s(midRamp) - startTime) ./ (endTime - startTime);
            contrastCurve(obj.time_s > endTime) = totalVolume_mm3(obj.time_s > endTime);
        end
    end
end
