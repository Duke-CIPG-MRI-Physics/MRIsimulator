classdef PhantomContext
    % PhantomContext
    %   Time base and shared waveform provider for the breast phantom
    %   components. Stores the reference time vector and precomputes shape
    %   parameters so downstream shapes can request full trajectories or
    %   indexed slices without duplicating arrays.

    properties (Access = private)
        time_s (1,:) double {mustBeFinite}
        heartA_mm (1,:) double {mustBeFinite, mustBeNonnegative}
        heartB_mm (1,:) double {mustBeFinite, mustBeNonnegative}
        heartC_mm (1,:) double {mustBeFinite, mustBeNonnegative}
        lungRadius_mm (1,:) double {mustBeFinite, mustBeNonnegative}
        lungHeight_mm (1,:) double {mustBeFinite, mustBeNonnegative}
        vesselContrast_mm3 (:,1) double {mustBeFinite, mustBeNonnegative}
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

        function [a_mm, b_mm, c_mm] = getHeartRadiiMm(obj, idx)
            % getHeartRadiiMm
            %   Accessor for the semi-axes of the beating heart ellipsoid.
            if nargin < 2 || isempty(idx)
                a_mm = obj.heartA_mm;
                b_mm = obj.heartB_mm;
                c_mm = obj.heartC_mm;
            else
                a_mm = obj.heartA_mm(idx);
                b_mm = obj.heartB_mm(idx);
                c_mm = obj.heartC_mm(idx);
            end
        end

        function r_mm = getLungRadiusMm(obj, idx)
            % getLungRadiusMm
            %   Accessor for lung radius over time.
            if nargin < 2 || isempty(idx)
                r_mm = obj.lungRadius_mm;
            else
                r_mm = obj.lungRadius_mm(idx);
            end
        end

        function h_mm = getLungHeightMm(obj, idx)
            % getLungHeightMm
            %   Accessor for lung height over time.
            if nargin < 2 || isempty(idx)
                h_mm = obj.lungHeight_mm;
            else
                h_mm = obj.lungHeight_mm(idx);
            end
        end

        function V_mm3 = getVesselContrastMm3(obj, idx)
            % getVesselContrastMm3
            %   Accessor for vessel contrast volume over time.
            if nargin < 2 || isempty(idx)
                V_mm3 = obj.vesselContrast_mm3;
            else
                V_mm3 = obj.vesselContrast_mm3(idx);
            end
        end
    end

    methods (Access = private)
        function [a_mm, b_mm, c_bm] = computeHeartGeometry(obj)
            heartOpts = struct('systFrac', 0.35, ...
                'q_ED', 50/27, ...
                'GLS_peak', -0.20, ...
                'GCS_peak', -0.25);
            numSamples = numel(obj.time_s);
            HR_bpm = 70 * ones(1, numSamples);
            EDV_ml = 150 * ones(1, numSamples);
            ESV_ml = 75  * ones(1, numSamples);

            [~, c_mm, a_mm, b_mm] = cardiac_ellipsoid_waveform(obj.time_s, HR_bpm, ...
                EDV_ml, ESV_ml, heartOpts);
        end

        function [R_mm, H_mm] = computeLungGeometry(obj)
            numSamples = numel(obj.time_s);
            f_bpm = 12 * ones(1, numSamples);
            VT_L = 0.4 * ones(1, numSamples);
            Vres_L = 0.8 * ones(1, numSamples);
            Vbase_L = 1.5 * ones(1, numSamples);
            bellyFrac = 0.6 * ones(1, numSamples);
            inspFrac = 0.4 * ones(1, numSamples);

            [~, R_mm, H_mm] = computeBreathingMotionEllipsoid(obj.time_s, f_bpm, VT_L, ...
                Vres_L, Vbase_L, bellyFrac, inspFrac);
        end

        function computeVesselContrast(obj)
            vesselRadius_mm = 2.5;
            total_vessel_length_mm = 100;

            totalVolume_mm3 = pi * vesselRadius_mm^2 * total_vessel_length_mm;
            V_contrast_mm3 = totalVolume_mm3 * ones(numel(obj.time_s), 1);

            obj.vesselContrast_mm3 = V_contrast_mm3;
        end
    end
end
