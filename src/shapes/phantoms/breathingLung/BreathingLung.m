classdef BreathingLung < SharedIntensityShapeGroup3D
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
    %       - Left and right lungs are created as AnalyticalEllipsoid3D
    %         components and supplied to the SharedIntensityShapeGroup3D
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
        R_mm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        H_mm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    methods
        function obj = BreathingLung(t_s, f_bpm, VT_L, Vres_L, Vbase_L, ...
                bellyFrac, inspFrac, lungSeparation_mm, intensity, center, rollPitchYaw)
            arguments
                t_s (1,:) double {mustBeReal, mustBeFinite}
                f_bpm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                VT_L (1,:) double {mustBeReal, mustBeFinite}
                Vres_L (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                Vbase_L (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
                bellyFrac (1,:) double {mustBeReal, mustBeFinite}
                inspFrac (1,:) double {mustBeReal, mustBeFinite}
                lungSeparation_mm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                intensity double {mustBeFinite} = 0.1
                center (1,3) double {mustBeFinite} = [0, 0, 0]
                rollPitchYaw (1,3) double {mustBeFinite} = [0 0 0]
            end

            [~, R_mm, H_mm] = computeBreathingMotionEllipsoid(t_s, f_bpm, VT_L, ...
                Vres_L, Vbase_L, bellyFrac, inspFrac);


            lungPosition_mm = R_mm + lungSeparation_mm;

            rightCenter = [lungPosition_mm(:), zeros(size(lungPosition_mm(:))), zeros(size(lungPosition_mm(:)))];
            leftCenter = [-lungPosition_mm(:), zeros(size(lungPosition_mm(:))), zeros(size(lungPosition_mm(:)))];


            rightLung = AnalyticalEllipsoid3D(R_mm, R_mm, H_mm, [], rightCenter, [0, 0, 0]);
            leftLung = AnalyticalEllipsoid3D(R_mm, R_mm, H_mm, [], leftCenter, [0, 0, 0]);

            rightLung.setTimeSamples(t_s);
            leftLung.setTimeSamples(t_s);

            obj@SharedIntensityShapeGroup3D([leftLung, rightLung], ...
                AnalyticalShape3D.empty, intensity, center, rollPitchYaw);

            obj.t_s = t_s;
            obj.f_bpm = f_bpm;
            obj.VT_L = VT_L;
            obj.Vres_L = Vres_L;
            obj.Vbase_L = Vbase_L;
            obj.bellyFrac = bellyFrac;
            obj.inspFrac = inspFrac;
            obj.R_mm = R_mm;
            obj.H_mm = H_mm;
        end

        function lungRadius = getLungRadiusMm(obj)
            lungRadius = obj.R_mm;
        end

        function lungHeight = getLungHeightMm(obj)
            lungHeight = obj.H_mm;
        end
    end
end
