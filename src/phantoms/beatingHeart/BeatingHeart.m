classdef BeatingHeart < AnalyticalEllipsoid3D
    % BeatingHeart
    %   Time-varying ellipsoidal heart model driven by cardiac_ellipsoid_waveform.
    %
    %   Constructor:
    %       heart = BeatingHeart(t_s, cardiacOpts, intensity, center, ...
    %           rollPitchYaw)
    %
    %   Inputs:
    %       t_s          : time samples [s]
    %       cardiacOpts  : struct with fields HR_bpm, EDV_ml, ESV_ml, and
    %                      optional parameters systFrac, q_ED, GLS_peak,
    %                      GCS_peak (see cardiac_ellipsoid_waveform)
    %       intensity    : (optional) shape intensity scaling
    %       center       : (optional) WORLD center [mm]
    %       rollPitchYaw : (optional) [roll pitch yaw] in deg
    %
    %   The constructor computes the semi-axes using cardiac_ellipsoid_waveform
    %   (which processes long signals in 15,625-sample chunks, ~0.001 Gb of
    %   temporaries) and initializes the parent AnalyticalEllipsoid3D with those
    %   axes.

    properties (Access = private)
        t_s (1,:) double {mustBeReal, mustBeFinite}
        HR_bpm double {mustBeReal, mustBeFinite, mustBeNonnegative}
        EDV_ml double {mustBeReal, mustBeFinite, mustBePositive}
        ESV_ml double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    methods
        function obj = BeatingHeart(t_s, cardiacOpts, intensity, center, rollPitchYaw)
            arguments
                t_s (1,:) double {mustBeReal, mustBeFinite}
                cardiacOpts struct
                intensity double {mustBeReal, mustBeFinite} = 1
                center double {mustBeReal, mustBeFinite} = [0 0 0]
                rollPitchYaw (1,3) double {mustBeReal, mustBeFinite} = [0 0 0]
            end

            heartParams = cardiac_ellipsoid_waveform(t_s, cardiacOpts);
            heartParams.pose = struct('center', struct('x_mm', center(:,1), ...
                'y_mm', center(:,2), ...
                'z_mm', center(:,3)), ...
                'roll_deg', rollPitchYaw(1), ...
                'pitch_deg', rollPitchYaw(2), ...
                'yaw_deg', rollPitchYaw(3));
            obj@AnalyticalEllipsoid3D(intensity, heartParams);

            obj.t_s = t_s;
            obj.HR_bpm = cardiacOpts.HR_bpm;
            obj.EDV_ml = cardiacOpts.EDV_ml;
            obj.ESV_ml = cardiacOpts.ESV_ml;
        end
    end
end
