classdef BeatingHeart < AnalyticalEllipsoid3D
    % BeatingHeart
    %   Time-varying ellipsoidal heart model driven by cardiac_ellipsoid_waveform.
    %
    %   Constructor:
    %       heart = BeatingHeart(t_s, HR_bpm, EDV_ml, ESV_ml, intensity, ...
    %           center, rollPitchYaw, opts)
    %
    %   Inputs:
    %       t_s          : time samples [s]
    %       HR_bpm       : heart rate waveform [beats/min]
    %       EDV_ml       : end-diastolic volume waveform [mL]
    %       ESV_ml       : end-systolic volume waveform [mL]
    %       intensity    : (optional) shape intensity scaling
    %       center       : (optional) WORLD center [mm]
    %       rollPitchYaw : (optional) [roll pitch yaw] in deg
    %       opts         : (optional) struct with fields systFrac, q_ED,
    %                      GLS_peak, GCS_peak (see cardiac_ellipsoid_waveform)
    %
    %   The constructor computes the semi-axes using cardiac_ellipsoid_waveform
    %   (which processes long signals in 15,625-sample chunks, ~0.001 Gb of
    %   temporaries) and initializes the parent AnalyticalEllipsoid3D with those
    %   axes.

    properties (Access = private)
        t_s (1,:) double {mustBeReal, mustBeFinite}
        HR_bpm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        EDV_ml (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
        ESV_ml (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    methods
        function obj = BeatingHeart(t_s, HR_bpm, EDV_ml, ESV_ml, intensity, ...
                center, rollPitchYaw, opts)
            arguments
                t_s (1,:) double {mustBeReal, mustBeFinite}
                HR_bpm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                EDV_ml (1,:) double {mustBeReal, mustBeFinite, mustBePositive}
                ESV_ml (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                intensity double {mustBeReal, mustBeFinite} = 1
                center double {mustBeReal, mustBeFinite} = [0 0 0]
                rollPitchYaw (1,3) double {mustBeReal, mustBeFinite} = [0 0 0]
                opts struct = struct()
            end

            [~, c_mm, a_mm, b_mm] = cardiac_ellipsoid_waveform(t_s, HR_bpm, ...
                EDV_ml, ESV_ml, opts);

            obj@AnalyticalEllipsoid3D(a_mm, b_mm, c_mm, intensity, center, ...
                rollPitchYaw);

            obj.t_s = t_s;
            obj.HR_bpm = HR_bpm;
            obj.EDV_ml = EDV_ml;
            obj.ESV_ml = ESV_ml;
        end
    end
end
