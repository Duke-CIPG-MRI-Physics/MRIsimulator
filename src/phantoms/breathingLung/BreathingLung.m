classdef BreathingLung < CompositeAnalyticalShape3D
    % BreathingLung
    %   Composite of two breathing lungs modeled as ellipsoids whose
    %   geometry follows lung_ellipsoid_waveform.
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
    %       lungSeparation_mm - additional spacing beyond lung radius [mm]
    %       intensity     - shape intensity (handled by parent)
    %       center        - composite center [mm]
    %       rollPitchYaw  - optional composite orientation [deg]
    %
    %   Notes:
    %       - Left and right lungs are created as AnalyticalEllipsoid3D
    %         components and supplied to the CompositeAnalyticalShape3D
    %         constructor as additive components.

    properties (Access = protected)
        t_s (1,:) double {mustBeReal, mustBeFinite}
        rightLung
        lefttLung
    end

    methods
        function obj = BreathingLung(t_s, pulmonaryOpts, lungSeparation_mm, intensity, shapeParameters)
            if nargin < 2 || isempty(shapeParameters)
                shapeParameters = BreathingLung.defaultLungParameters();
            end

            % Calculate breathing waveform
            ellipsoidParams = lung_ellipsoid_waveform(t_s, pulmonaryOpts);

            lungParams = struct('a_mm', ellipsoidParams.R_mm, ...
                'b_mm', ellipsoidParams.R_mm, ...
                'c_mm', ellipsoidParams.H_mm);
            
            rightCenter_mm = [ellipsoidParams.lungPosition_mm(:), ...
                zeros(numel(ellipsoidParams.lungPosition_mm), 1), ...
                zeros(numel(ellipsoidParams.lungPosition_mm), 1)];
            leftCenter_mm = [-ellipsoidParams.lungPosition_mm(:), ...
                zeros(numel(ellipsoidParams.lungPosition_mm), 1), ...
                zeros(numel(ellipsoidParams.lungPosition_mm), 1)];
            rightParams = lungParams;
            rightParams.pose = struct('center', struct('x_mm', rightCenter_mm(:,1), ...
                'y_mm', rightCenter_mm(:,2), ...
                'z_mm', rightCenter_mm(:,3)), ...
                'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);
            leftParams = lungParams;
            leftParams.pose = struct('center', struct('x_mm', leftCenter_mm(:,1), ...
                'y_mm', leftCenter_mm(:,2), ...
                'z_mm', leftCenter_mm(:,3)), ...
                'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);

            rightLung = AnalyticalEllipsoid3D([], rightParams);
            leftLung = AnalyticalEllipsoid3D([], leftParams);

            obj@CompositeAnalyticalShape3D([leftLung, rightLung], ...
                AnalyticalShape3D.empty, intensity, shapeParameters);

            obj.t_s = t_s;
            obj.rightLung = rightLung;
            obj.leftLung = leftLung;
        end

        function leftLung = getLeftLung(obj)
            leftLung = obj.leftLung;
        end

        function rightLung = getRightLung(obj)
            rightLung = obj.rightLung;
        end
    end
end
