classdef PhantomContext
    % PhantomContext
    %   Time base and shared waveform provider for the breast phantom
    %   components. Downstream shapes request geometry through function
    %   handles so parameters are computed on demand rather than stored in
    %   memory-heavy arrays.

    properties (Access = private)
        time_s (1,:) double {mustBeFinite}
    end

    methods
        function obj = PhantomContext(t_s)
            % PhantomContext
            %   Provider for time-varying parameters used by breast phantom
            %   components.
            %
            %   obj = PhantomContext(t_s)
            %       t_s : column vector of time samples [s]
            arguments
                t_s (:,1) double {mustBeFinite}
            end

            obj.time_s = t_s(:).';
        end

        function t = getTime(obj, idx)
            % getTime
            %   Return the full time vector or a slice specified by idx.
            if nargin < 2 || isempty(idx)
                t = obj.time_s;
            else
                t = obj.time_s(idx);
            end
        end

        function [aFcn, bFcn, cFcn] = getHeartRadiiMm(obj)
            % getHeartRadiiMm
            %   Return function handles that compute the semi-axes of the
            %   beating heart ellipsoid for a requested time vector. Handles
            %   accept a row vector of time samples and return matching
            %   trajectories.

            aFcn = @(t) obj.computeHeartComponent(t, 'a');
            bFcn = @(t) obj.computeHeartComponent(t, 'b');
            cFcn = @(t) obj.computeHeartComponent(t, 'c');
        end

        function radiusFcn = getLungRadiusMm(obj)
            % getLungRadiusMm
            %   Return a function handle that computes lung radius for a
            %   given time vector.

            radiusFcn = @(t) obj.computeLungGeometry(obj.normalizeTime(t), 'radius');
        end

        function heightFcn = getLungHeightMm(obj)
            % getLungHeightMm
            %   Return a function handle that computes lung height for a
            %   given time vector.

            heightFcn = @(t) obj.computeLungGeometry(obj.normalizeTime(t), 'height');
        end

        function contrastFcn = getVesselContrastMm3(obj)
            % getVesselContrastMm3
            %   Return a function handle that computes vessel contrast
            %   volume for a given time vector.

            contrastFcn = @(t) obj.computeVesselContrast(obj.normalizeTime(t));
        end
    end

    methods (Access = private)
        function tRow = normalizeTime(obj, t_s)
            if nargin < 2 || isempty(t_s)
                tRow = obj.time_s;
            else
                validateattributes(t_s, {'double'}, {'finite'});
                tRow = t_s(:).';
            end
        end

        function [a_mm, b_mm, c_mm] = computeHeartGeometry(obj, t_s)
            heartOpts = struct('systFrac', 0.35, ...
                'q_ED', 50/27, ...
                'GLS_peak', -0.20, ...
                'GCS_peak', -0.25);
            numSamples = numel(t_s);
            HR_bpm = 70 * ones(1, numSamples);
            EDV_ml = 150 * ones(1, numSamples);
            ESV_ml = 75  * ones(1, numSamples);

            [~, c_mm, a_mm, b_mm] = cardiac_ellipsoid_waveform(t_s, HR_bpm, ...
                EDV_ml, ESV_ml, heartOpts);
        end

        function component = computeHeartComponent(obj, t_s, whichComponent)
            [a_mm, b_mm, c_mm] = obj.computeHeartGeometry(obj.normalizeTime(t_s));

            switch whichComponent
                case 'a'
                    component = a_mm;
                case 'b'
                    component = b_mm;
                case 'c'
                    component = c_mm;
                otherwise
                    error('PhantomContext:UnknownHeartComponent', ...
                        'Unknown heart component: %s', whichComponent);
            end
        end

        function geometry = computeLungGeometry(obj, t_s, whichComponent)
            numSamples = numel(t_s);
            f_bpm = 12 * ones(1, numSamples);
            VT_L = 0.4 * ones(1, numSamples);
            Vres_L = 0.8 * ones(1, numSamples);
            Vbase_L = 1.5 * ones(1, numSamples);
            bellyFrac = 0.6 * ones(1, numSamples);
            inspFrac = 0.4 * ones(1, numSamples);

            [~, R_mm, H_mm] = computeBreathingMotionEllipsoid(t_s, f_bpm, VT_L, ...
                Vres_L, Vbase_L, bellyFrac, inspFrac);

            switch whichComponent
                case 'radius'
                    geometry = R_mm;
                case 'height'
                    geometry = H_mm;
                otherwise
                    error('PhantomContext:UnknownLungComponent', ...
                        'Unknown lung component: %s', whichComponent);
            end
        end

        function V_contrast_mm3 = computeVesselContrast(obj, t_s)
            vesselRadius_mm = 2.5;
            total_vessel_length_mm = 100;

            totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
            V_contrast_mm3 = totalVolume_mm3 * ones(numel(t_s), 1);
        end
    end
end
