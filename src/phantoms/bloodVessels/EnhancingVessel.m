classdef EnhancingVessel < MultipleMaterialPhantom
    % EnhancingVessel
    %   Vessel phantom built from two AnalyticalCylinder3D objects: an
    %   enhancing segment (with contrast) and an unenhanced segment. The
    %   total length remains fixed while the enhanced segment grows or
    %   shrinks according to a supplied contrast volume curve. The object
    %   behaves as a MultipleMaterialPhantom so it can be added directly to
    %   other phantoms.

    properties (Access = private)
        enhancedVessel  % No type constraint
        unenhancedVessel  % No type constraint
        totalLength_mm double {mustBePositive} = 1;
        enhancingVesselCenter (1,3) double = [0 0 0];
        vesselRadius_mm double {mustBePositive} = 1;
        contrastVolume_mm3 (:,1) double {mustBeFinite, mustBeNonnegative}

        time_s (:,1) double {mustBeFinite}
        rollPitchYaw (1,3) double = [0 0 0];
    end

    methods
        function obj = EnhancingVessel(t_s, totalLength_mm, ...
                enhancedIntensity, unenhancedIntensity, ...
                vesselRadius_mm, V_contrast_mm3, vesselCenter, rollPitchYaw)
            % Constructor
            %   obj = EnhancingVessel(t_s, totalLength_mm, enhancedIntensity, unenhancedIntensity, vesselRadius_mm, V_contrast_mm3)
            %   obj = EnhancingVessel(..., centerlineCenter)
            %   obj = EnhancingVessel(..., centerlineCenter, rollPitchYaw)
            %
            %   t_s             : time samples [s]
            %   totalLength_mm:  scalar or vector total vessel length [mm]
            %   enhancedIntensity:  intensity assigned to the contrast-enhanced segment
            %   unenhancedIntensity:intensity assigned to the baseline segment
            %   vesselRadius_mm:    radius over time (scalar or vector matching t_s)
            %   V_contrast_mm3 :    contrast volume over time [mm^3]
            %   centerlineCenter:   [1x3] WORLD coords for the combined vessel center
            %   rollPitchYaw:       [1x3] WORLD orientation of vessel body axes (deg)

            arguments
                t_s (:,1) double {mustBeFinite}
                totalLength_mm double {mustBePositive}
                enhancedIntensity double
                unenhancedIntensity double
                vesselRadius_mm double {mustBeNonnegative}
                V_contrast_mm3 (:,1) double {mustBeFinite, mustBeNonnegative}
                vesselCenter (1,3) double = [0 0 0];
                rollPitchYaw (1,3) double = [0 0 0];
            end

            obj@MultipleMaterialPhantom();

            obj.time_s = t_s;
            obj.totalLength_mm = obj.ensureTotalLengthVector(totalLength_mm, numel(t_s));
            obj.vesselRadius_mm = obj.ensureRadiusVector(vesselRadius_mm, numel(t_s));
            obj.contrastVolume_mm3 = obj.ensureVectorSize(V_contrast_mm3, numel(t_s));
            obj.rollPitchYaw = rollPitchYaw;

            [enhancedLength_mm, unenhancedLength_mm] = obj.calcVesselLengths(obj.time_s, obj.contrastVolume_mm3, obj.vesselRadius_mm);
            [enhancedCenter_mm, unenhancedCenter_mm] = obj.calcVesselCenter(vesselCenter, ...
                enhancedLength_mm, unenhancedLength_mm);

            enhancedParams = struct('radius_mm', obj.vesselRadius_mm, ...
                'length_mm', enhancedLength_mm);
            enhancedParams.pose = struct('center', struct('x_mm', enhancedCenter_mm(:,1), ...
                'y_mm', enhancedCenter_mm(:,2), ...
                'z_mm', enhancedCenter_mm(:,3)), ...
                'roll_deg', rollPitchYaw(1), ...
                'pitch_deg', rollPitchYaw(2), ...
                'yaw_deg', rollPitchYaw(3));
            obj.enhancedVessel = AnalyticalCylinder3D(enhancedIntensity, enhancedParams);

            unenhancedParams = struct('radius_mm', obj.vesselRadius_mm, ...
                'length_mm', unenhancedLength_mm);
            unenhancedParams.pose = struct('center', struct('x_mm', unenhancedCenter_mm(:,1), ...
                'y_mm', unenhancedCenter_mm(:,2), ...
                'z_mm', unenhancedCenter_mm(:,3)), ...
                'roll_deg', rollPitchYaw(1), ...
                'pitch_deg', rollPitchYaw(2), ...
                'yaw_deg', rollPitchYaw(3));
            obj.unenhancedVessel = AnalyticalCylinder3D(unenhancedIntensity, unenhancedParams);

            obj.setShapes([obj.unenhancedVessel, obj.enhancedVessel]);
        end

        function updateTimeArray(obj, t_s, contrastVolume_mm3, vesselRadius_mm)
            % updateTimeArray
            %   Recompute the enhanced/unenhanced trajectories for a new time
            %   vector and update the paired cylindrical vessels.
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
                contrastVolume_mm3 (:,1) double {mustBeFinite, mustBeNonnegative} = obj.contrastVolume_mm3
                vesselRadius_mm double {mustBePositive} = obj.vesselRadius_mm
            end

            obj.time_s = t_s;
            obj.totalLength_mm = obj.ensureTotalLengthVector(obj.totalLength_mm, numel(t_s));
            obj.contrastVolume_mm3 = obj.ensureVectorSize(contrastVolume_mm3, numel(t_s), ...
                'EnhancingVessel:ContrastSizeMismatch', ...
                'V_contrast_mm3 must match the length of t_s.');
            obj.vesselRadius_mm = obj.ensureRadiusVector(vesselRadius_mm, numel(t_s));

            obj.updateVesselLengths(obj.time_s, obj.contrastVolume_mm3, obj.vesselRadius_mm);
            obj.updateVesselCenter(obj.enhancingVesselCenter);
        end

        function setCenterlineCenter(obj, newCenter)
            % setCenterlineCenter
            %   Move the combined vessel while preserving the relative
            %   placement of the enhanced/unenhanced segments.
            arguments
                obj
                newCenter (1,3) double
            end
            obj.enhancingVesselCenter = double(newCenter(:)).';
            obj.updateVesselCenter(obj.enhancingVesselCenter);
        end

        function [enhancedVessel, unenhancedVessel] = getVessels(obj)
            % getVessels
            %   Return handles to the enhanced and unenhanced cylinder
            %   objects for downstream manipulation or inspection.
            enhancedVessel = obj.enhancedVessel;
            unenhancedVessel = obj.unenhancedVessel;
        end
    end

    methods (Access = private)
        function axisUnit = bodyZAxisInWorld(obj, rpy)
            % bodyZAxisInWorld
            %   Convenience helper that returns the vessel body +z axis in
            %   world coordinates. If an explicit roll/pitch/yaw is not
            %   supplied, the stored vessel orientation is used; if the
            %   vessels have not been constructed yet, a zero-orientation is
            %   assumed.
            if nargin < 2 || isempty(rpy)
                if ~isempty(obj.enhancedVessel) && isvalid(obj.enhancedVessel)
                    params = obj.enhancedVessel.getShapeParameters();
                    pose = obj.extractPose(params, 0);
                    rpy = [pose.roll_deg(1), pose.pitch_deg(1), pose.yaw_deg(1)];
                else
                    rpy = obj.rollPitchYaw;
                end
            end
            R = obj.rotationMatrixFromRPY(rpy);
            axisUnit = (R * [0; 0; 1]).';
            axisUnit = axisUnit ./ norm(axisUnit);
        end

        function [enhancedLength_mm, unenhancedLength_mm] = calcVesselLengths(obj, t_s, V_contrast_mm3, vesselRadius_mm)
            % calcVesselLengths
            %   Compute the time-varying lengths of the enhanced and
            %   unenhanced segments given a contrast volume curve and vessel
            %   radius.
            enhancedLength_mm = computeContrastWashIn(t_s, vesselRadius_mm, V_contrast_mm3);
            enhancedLength_mm = min(enhancedLength_mm, obj.totalLength_mm);

            unenhancedLength_mm = obj.totalLength_mm - enhancedLength_mm;
            unenhancedLength_mm = max(unenhancedLength_mm, 0);
        end

        function updateVesselLengths(obj, t_s, V_contrast_mm3, vesselRadius_mm)
            % updateVesselLengths
            %   Recompute segment lengths and apply them to the vessel
            %   shapes.
            [enhancedLength_mm, unenhancedLength_mm] = obj.calcVesselLengths(t_s, V_contrast_mm3, vesselRadius_mm);
            enhancedParams = obj.enhancedVessel.getShapeParameters();
            enhancedParams.radius_mm = vesselRadius_mm;
            enhancedParams.length_mm = enhancedLength_mm;
            obj.enhancedVessel.setShapeParameters(enhancedParams);

            unenhancedParams = obj.unenhancedVessel.getShapeParameters();
            unenhancedParams.radius_mm = vesselRadius_mm;
            unenhancedParams.length_mm = unenhancedLength_mm;
            obj.unenhancedVessel.setShapeParameters(unenhancedParams);
        end

        function [enhancedCenter, unenhancedCenter] = calcVesselCenter(obj, vesselCenter, enhancedLength_mm, unenhancedLength_mm)
            % calcVesselCenter
            %   Calculate the centers of the enhanced and unenhanced
            %   segments 
            axisUnit = obj.bodyZAxisInWorld(obj.rollPitchYaw);

            enhancedCenter = vesselCenter - 0.5 .* unenhancedLength_mm .* axisUnit;
            unenhancedCenter = vesselCenter + 0.5 .* enhancedLength_mm .* axisUnit;
        end

        function updateVesselCenters(obj, vesselCenter)
            % updateVesselCenter
            %   Reposition the enhanced and unenhanced vessels while
            %   maintaining their relative spacing along the body z-axis.
            enhancedParams = obj.enhancedVessel.getShapeParameters();
            unenhancedParams = obj.unenhancedVessel.getShapeParameters();

            enhancedLength_mm = enhancedParams.length_mm;
            unenhancedLength_mm = unenhancedParams.length_mm;

            [enhancedCenter, unenhancedCenter] = obj.calcVesselCenter(vesselCenter, ...
                enhancedLength_mm, unenhancedLength_mm);

            enhancedParams.pose.center.x_mm = enhancedCenter(:,1);
            enhancedParams.pose.center.y_mm = enhancedCenter(:,2);
            enhancedParams.pose.center.z_mm = enhancedCenter(:,3);
            obj.enhancedVessel.setShapeParameters(enhancedParams);

            unenhancedParams.pose.center.x_mm = unenhancedCenter(:,1);
            unenhancedParams.pose.center.y_mm = unenhancedCenter(:,2);
            unenhancedParams.pose.center.z_mm = unenhancedCenter(:,3);
            obj.unenhancedVessel.setShapeParameters(unenhancedParams);

            phantomPose = obj.getShapeParameters();
            phantomPose.pose.center.x_mm = vesselCenter(:,1);
            phantomPose.pose.center.y_mm = vesselCenter(:,2);
            phantomPose.pose.center.z_mm = vesselCenter(:,3);
            obj.setShapeParameters(phantomPose);
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

        function lengthVec = ensureTotalLengthVector(~, totalLengthInput, targetLength)
            if isscalar(totalLengthInput)
                lengthVec = totalLengthInput * ones(targetLength, 1);
            else
                if numel(totalLengthInput) ~= targetLength
                    error('EnhancingVessel:TotalLengthSizeMismatch', ...
                        'totalLength_mm must be scalar or match length of t_s.');
                end
                lengthVec = totalLengthInput(:);
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
