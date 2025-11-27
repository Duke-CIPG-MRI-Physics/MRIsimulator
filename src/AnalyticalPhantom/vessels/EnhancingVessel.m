classdef EnhancingVessel < MultipleMaterialPhantom
    % EnhancingVessel
    %   Vessel phantom built from two AnalyticalEllipticalCylinder3D objects:
    %   an enhancing segment (with contrast) and an unenhanced segment. The
    %   total length remains fixed while the enhanced segment grows or
    %   shrinks according to a supplied contrast volume curve. The object
    %   behaves as a MultipleMaterialPhantom so it can be added directly to
    %   other phantoms.

    properties (Access = private)
        enhancedVessel (1,1) AnalyticalEllipticalCylinder3D
        unenhancedVessel (1,1) AnalyticalEllipticalCylinder3D
        totalLength_mm (1,1) double {mustBePositive}
        centerlineCenter (1,3) double = [0 0 0];
        vesselRadius_mm double {mustBePositive}
        radiusInput
        contrastVolume_mm3 (:,1) double {mustBeFinite, mustBeNonnegative}

        time_s (:,1) double {mustBeFinite}
        enhancedLengths_mm (:,1) double {mustBeFinite, mustBeNonnegative}
        unenhancedLengths_mm (:,1) double {mustBeFinite, mustBeNonnegative}
        enhancedCenters_mm (:,3) double {mustBeFinite}
        unenhancedCenters_mm (:,3) double {mustBeFinite}
    end

    methods
        function obj = EnhancingVessel(t_s, totalLength_mm, enhancedIntensity, unenhancedIntensity, vesselRadius_mm, V_contrast_mm3, centerlineCenter, rollPitchYaw)
            % Constructor
            %   obj = EnhancingVessel(t_s, totalLength_mm, enhancedIntensity, unenhancedIntensity, vesselRadius_mm, V_contrast_mm3)
            %   obj = EnhancingVessel(..., centerlineCenter)
            %   obj = EnhancingVessel(..., centerlineCenter, rollPitchYaw)
            %
            %   t_s             : time samples [s]
            %   totalLength_mm:     scalar total vessel length [mm]
            %   enhancedIntensity:  intensity assigned to the contrast-enhanced segment
            %   unenhancedIntensity:intensity assigned to the baseline segment
            %   vesselRadius_mm:    radius over time (scalar or vector matching t_s)
            %   V_contrast_mm3 :    contrast volume over time [mm^3]
            %   centerlineCenter:   [1x3] WORLD coords for the combined vessel center
            %   rollPitchYaw:       [1x3] WORLD orientation of vessel body axes (deg)

            arguments
                t_s (:,1) double {mustBeFinite}
                totalLength_mm (1,1) double {mustBePositive}
                enhancedIntensity double
                unenhancedIntensity double
                vesselRadius_mm double {mustBePositive}
                V_contrast_mm3 (:,1) double {mustBeFinite, mustBeNonnegative}
                centerlineCenter (1,3) double = [0 0 0];
                rollPitchYaw (1,3) double = [0 0 0];
            end

            obj@MultipleMaterialPhantom();

            obj.totalLength_mm = totalLength_mm;
            obj.centerlineCenter = centerlineCenter;
            obj.time_s = t_s;
            obj.radiusInput = vesselRadius_mm;
            obj.contrastVolume_mm3 = V_contrast_mm3(:);

            obj.updateGeometrySeries(rollPitchYaw);

            obj.enhancedVessel = AnalyticalEllipticalCylinder3D(obj.vesselRadius_mm, obj.vesselRadius_mm, ...
                obj.enhancedLengths_mm, enhancedIntensity, obj.enhancedCenters_mm, rollPitchYaw);
            obj.unenhancedVessel = AnalyticalEllipticalCylinder3D(obj.vesselRadius_mm, obj.vesselRadius_mm, ...
                obj.unenhancedLengths_mm, unenhancedIntensity, obj.unenhancedCenters_mm, rollPitchYaw);

            obj.setShapes([obj.unenhancedVessel, obj.enhancedVessel]);
        end

        function updateTimeArray(obj, t_s, V_contrast_mm3, vesselRadius_mm)
            % updateTimeArray
            %   Recompute the enhanced/unenhanced trajectories for a new time
            %   vector and update the paired elliptical cylinders.
            %
            %   Inputs:
            %       t_s             : time samples [s]
            %       V_contrast_mm3  : contrast volume over time [mm^3]
            %                         (optional; defaults to stored values)
            %       vesselRadius_mm : radius over time [mm]
            %                         (optional; defaults to stored values)

            arguments
                obj
                t_s (:,1) double {mustBeFinite}
                V_contrast_mm3 (:,1) double {mustBeFinite, mustBeNonnegative} = obj.contrastVolume_mm3
                vesselRadius_mm double {mustBePositive} = obj.radiusInput
            end

            obj.time_s = t_s;
            obj.contrastVolume_mm3 = obj.ensureVectorSize(V_contrast_mm3, numel(t_s), ...
                'EnhancingVessel:ContrastSizeMismatch', ...
                'V_contrast_mm3 must match the length of t_s.');
            obj.radiusInput = vesselRadius_mm;

            obj.updateGeometrySeries();

            obj.enhancedVessel.setLength(obj.enhancedLengths_mm);
            obj.unenhancedVessel.setLength(obj.unenhancedLengths_mm);

            obj.enhancedVessel.setA(obj.vesselRadius_mm);
            obj.enhancedVessel.setB(obj.vesselRadius_mm);
            obj.unenhancedVessel.setA(obj.vesselRadius_mm);
            obj.unenhancedVessel.setB(obj.vesselRadius_mm);

            obj.enhancedVessel.setCenter(obj.enhancedCenters_mm);
            obj.unenhancedVessel.setCenter(obj.unenhancedCenters_mm);
        end

        function setCenterlineCenter(obj, newCenter)
            % setCenterlineCenter
            %   Move the combined vessel while preserving the relative
            %   placement of the enhanced/unenhanced segments.
            arguments
                obj
                newCenter (1,3) double
            end
            obj.centerlineCenter = double(newCenter(:)).';
            obj.refreshSegmentCenters();
            obj.enhancedVessel.setCenter(obj.enhancedCenters_mm);
            obj.unenhancedVessel.setCenter(obj.unenhancedCenters_mm);
        end

        function [enhancedVessel, unenhancedVessel] = getVessels(obj)
            enhancedVessel = obj.enhancedVessel;
            unenhancedVessel = obj.unenhancedVessel;
        end
    end

    methods (Access = private)
        function axisUnit = bodyZAxisInWorld(obj, rpy)
            if nargin < 2 || isempty(rpy)
                rpy = obj.enhancedVessel.getRollPitchYaw();
            end
            R = obj.rotationMatrixFromRPY(rpy);
            axisUnit = (R * [0; 0; 1]).';
            axisUnit = axisUnit ./ norm(axisUnit);
        end

        function [enhancedCenters, unenhancedCenters] = computeSegmentCenters(obj, axisUnit)
            enhancedCenters = obj.centerlineCenter - 0.5 * obj.unenhancedLengths_mm .* axisUnit;
            unenhancedCenters = obj.centerlineCenter + 0.5 * obj.enhancedLengths_mm .* axisUnit;
        end

        function refreshSegmentCenters(obj)
            axisUnit = obj.bodyZAxisInWorld();
            [obj.enhancedCenters_mm, obj.unenhancedCenters_mm] = obj.computeSegmentCenters(axisUnit);
        end

        function updateGeometrySeries(obj, rollPitchYaw)
            if nargin < 2
                rollPitchYaw = [];
            end

            obj.vesselRadius_mm = obj.ensureRadiusVector(obj.radiusInput, numel(obj.time_s));

            obj.enhancedLengths_mm = computeContrastWashIn(obj.time_s, obj.vesselRadius_mm, obj.contrastVolume_mm3);
            obj.enhancedLengths_mm = min(obj.enhancedLengths_mm, obj.totalLength_mm);
            obj.unenhancedLengths_mm = max(obj.totalLength_mm - obj.enhancedLengths_mm, eps);
            obj.enhancedLengths_mm = max(obj.enhancedLengths_mm, eps);

            axisUnit = obj.bodyZAxisInWorld(rollPitchYaw);
            [obj.enhancedCenters_mm, obj.unenhancedCenters_mm] = obj.computeSegmentCenters(axisUnit);
        end

        function radiusVec = ensureRadiusVector(~, radiusInput, targetLength)
            if isscalar(radiusInput)
                radiusVec = radiusInput * ones(targetLength, 1);
            else
                if numel(radiusInput) ~= targetLength
                    error('EnhancingVessel:RadiusSizeMismatch', ...
                        'vesselRadius_mm must be scalar or match length of t_s.');
                end
                radiusVec = radiusInput(:);
            end
        end

        function vec = ensureVectorSize(~, inputVec, targetLength, id, msg)
            inputVec = inputVec(:);
            if numel(inputVec) ~= targetLength
                error(id, msg);
            end
            vec = inputVec;
        end

        function R = rotationMatrixFromRPY(~, rpy)
            r = deg2rad(rpy(1));
            p = deg2rad(rpy(2));
            y = deg2rad(rpy(3));

            Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
            Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
            Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];

            R = Rz * Ry * Rx;
        end

    end
end
