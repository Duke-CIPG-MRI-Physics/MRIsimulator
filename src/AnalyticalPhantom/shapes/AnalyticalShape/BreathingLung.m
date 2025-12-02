classdef BreathingLung < CompositeAnalyticalShape3D
    % BreathingLung
    %   Composite of two breathing lungs modeled as ellipsoids whose
    %   geometry follows computeBreathingMotionEllipsoid.
    %
    %   Constructor:
    %       obj = BreathingLung(t_s, f_bpm, VT_L, Vres_L, Vbase_L, ...
    %           bellyFrac, inspFrac, maxHeartDim_mm, intensity, center, rollPitchYaw)
    %
    %   Inputs:
    %       t_s           - time [s]
    %       f_bpm         - breathing frequency [breaths/min]
    %       VT_L          - tidal volume [L]
    %       Vres_L        - residual volume [L]
    %       Vbase_L       - baseline lung volume [L]
    %       bellyFrac     - belly-breathing fraction [0,1]
    %       inspFrac      - inspiratory fraction of cycle (0,1)
    %       maxHeartDim_mm- maximal heart dimension [mm] for lung spacing
    %       intensity     - shape intensity (handled by parent)
    %       center        - composite center [mm]
    %       rollPitchYaw  - optional composite orientation [deg]
    %
    %   Notes:
    %       - Left and right lungs are stored as private AnalyticalEllipsoid3D
    %         components and supplied to the CompositeAnalyticalShape3D
    %         constructor as additive components.

    properties (Access = private)
        t_s (1,:) double {mustBeReal, mustBeFinite}
        f_bpm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        VT_L (1,:) double {mustBeReal, mustBeFinite}
        Vres_L (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        Vbase_L (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        bellyFrac (1,:) double {mustBeReal, mustBeFinite}
        inspFrac (1,:) double {mustBeReal, mustBeFinite}
        maxHeartDim_mm (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        leftLung (1,1) AnalyticalEllipsoid3D
        rightLung (1,1) AnalyticalEllipsoid3D
        lungRadius_mm (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        lungHeight_mm (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        lungSeparation_mm (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        maxLungSize_mm (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    methods
        function obj = BreathingLung(t_s, f_bpm, VT_L, Vres_L, Vbase_L, ...
                bellyFrac, inspFrac, maxHeartDim_mm, intensity, center, rollPitchYaw)
            arguments
                t_s (1,:) double {mustBeReal, mustBeFinite}
                f_bpm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                VT_L (1,:) double {mustBeReal, mustBeFinite}
                Vres_L (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                Vbase_L (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
                bellyFrac (1,:) double {mustBeReal, mustBeFinite}
                inspFrac (1,:) double {mustBeReal, mustBeFinite}
                maxHeartDim_mm (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                intensity double {mustBeFinite} = 0.1
                center (1,3) double {mustBeFinite} = [0, 0, 0]
                rollPitchYaw (1,3) double {mustBeFinite} = [0 0 0]
            end

            obj.t_s = t_s;
            obj.f_bpm = f_bpm;
            obj.VT_L = VT_L;
            obj.Vres_L = Vres_L;
            obj.Vbase_L = Vbase_L;
            obj.bellyFrac = bellyFrac;
            obj.inspFrac = inspFrac;
            obj.maxHeartDim_mm = maxHeartDim_mm;

            [~, R_m, H_m] = computeBreathingMotionEllipsoid(t_s, f_bpm, VT_L, ...
                Vres_L, Vbase_L, bellyFrac, inspFrac);

            R_mm = 1000 .* R_m;
            H_mm = 1000 .* H_m;

            obj.lungRadius_mm = max(R_mm(:));
            obj.lungHeight_mm = max(H_mm(:));
            obj.maxLungSize_mm = max(obj.lungRadius_mm, obj.lungHeight_mm);
            obj.lungSeparation_mm = obj.lungRadius_mm + maxHeartDim_mm + 2;

            rightCenter = [obj.lungSeparation_mm, 0, 0];
            leftCenter = [-obj.lungSeparation_mm, 0, 0];

            obj.rightLung = AnalyticalEllipsoid3D(R_mm, R_mm, H_mm, [], rightCenter, [0, 95, 0]);
            obj.leftLung = AnalyticalEllipsoid3D(R_mm, R_mm, H_mm, [], leftCenter, [0, 85, 0]);

            obj@CompositeAnalyticalShape3D([obj.leftLung, obj.rightLung], ...
                AnalyticalShape3D.empty, intensity, center, rollPitchYaw);
        end

        function lungRadius = getLungRadiusMm(obj)
            lungRadius = obj.lungRadius_mm;
        end

        function lungHeight = getLungHeightMm(obj)
            lungHeight = obj.lungHeight_mm;
        end

        function separation = getLungSeparationMm(obj)
            separation = obj.lungSeparation_mm;
        end

        function maxSize = getMaxLungSizeMm(obj)
            maxSize = obj.maxLungSize_mm;
        end
    end
end
