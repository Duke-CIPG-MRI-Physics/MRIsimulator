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
    %   and initializes the parent AnalyticalEllipsoid3D with those axes.

    properties (Access = private)
        t_s (1,:) double {mustBeReal, mustBeFinite}
        HR_bpm
        EDV_ml
        ESV_ml
        waveformOpts struct = struct()
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

            obj.waveformOpts = opts;
            obj.HR_bpm = obj.normalizeGeometryInput(HR_bpm, ...
                @(v) validateattributes(v, {'double'}, {'real', 'finite', 'nonnegative'}), ...
                true, 'HR_bpm');
            obj.EDV_ml = obj.normalizeGeometryInput(EDV_ml, ...
                @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'}), ...
                true, 'EDV_ml');
            obj.ESV_ml = obj.normalizeGeometryInput(ESV_ml, ...
                @(v) validateattributes(v, {'double'}, {'real', 'finite', 'nonnegative'}), ...
                true, 'ESV_ml');

            obj@AnalyticalEllipsoid3D([], [], [], intensity, center, rollPitchYaw);

            obj.t_s = t_s;
            obj.setTimeSamples(t_s);
            obj.updateGeometryFromWaveform();
        end

        function setTimeSamples(obj, t_s)
            % setTimeSamples  Override to keep stored time base in sync and
            % refresh waveform-driven geometry.
            arguments
                obj
                t_s (1,:) double {mustBeReal, mustBeFinite}
            end

            obj.t_s = t_s;
            setTimeSamples@AnalyticalEllipsoid3D(obj, t_s);
            obj.updateGeometryFromWaveform();
        end

        function setEDVml(obj, newEDV_ml, opts)
            % setEDVml  Update end-diastolic volume (numeric or @(t)->EDV_ml).
            %   Function handles should accept a 1xN time vector and return a
            %   scalar or vector of matching length. Use Cache=false to force
            %   re-evaluation on every geometry query.
            arguments
                obj
                newEDV_ml
                opts.Cache logical = true
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'positive'});
            edvSpec = obj.normalizeGeometryInput(newEDV_ml, validator, opts.Cache, 'EDV_ml');

            if ~isequal(obj.EDV_ml, edvSpec)
                obj.EDV_ml = edvSpec;
                obj.updateGeometryFromWaveform();
            end
        end

        function setESVml(obj, newESV_ml, opts)
            % setESVml  Update end-systolic volume (numeric or @(t)->ESV_ml).
            %   Function handles should accept a 1xN time vector and return a
            %   scalar or vector of matching length.
            arguments
                obj
                newESV_ml
                opts.Cache logical = true
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'nonnegative'});
            esvSpec = obj.normalizeGeometryInput(newESV_ml, validator, opts.Cache, 'ESV_ml');

            if ~isequal(obj.ESV_ml, esvSpec)
                obj.ESV_ml = esvSpec;
                obj.updateGeometryFromWaveform();
            end
        end

        function setHeartRateBpm(obj, newHR_bpm, opts)
            % setHeartRateBpm  Update HR waveform (numeric or @(t)->HR_bpm).
            arguments
                obj
                newHR_bpm
                opts.Cache logical = true
            end

            validator = @(v) validateattributes(v, {'double'}, {'real', 'finite', 'nonnegative'});
            hrSpec = obj.normalizeGeometryInput(newHR_bpm, validator, opts.Cache, 'HR_bpm');

            if ~isequal(obj.HR_bpm, hrSpec)
                obj.HR_bpm = hrSpec;
                obj.updateGeometryFromWaveform();
            end
        end
    end

    methods (Access = private)
        function updateGeometryFromWaveform(obj)
            edv_ml = obj.requireScalarOrSize(obj.EDV_ml, obj.t_s, 'EDV_ml');
            esv_ml = obj.requireScalarOrSize(obj.ESV_ml, obj.t_s, 'ESV_ml');
            hr_bpm = obj.requireScalarOrSize(obj.HR_bpm, obj.t_s, 'HR_bpm');

            [~, c_mm, a_mm, b_mm] = cardiac_ellipsoid_waveform(obj.t_s, hr_bpm, ...
                edv_ml, esv_ml, obj.waveformOpts);

            obj.setA(a_mm);
            obj.setB(b_mm);
            obj.setC(c_mm);
        end
    end
end
