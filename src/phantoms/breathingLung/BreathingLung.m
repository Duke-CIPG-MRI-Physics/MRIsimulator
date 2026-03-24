classdef BreathingLung < CompositeAnalyticalShape3D
    % BreathingLung
    %   Composite of two breathing lungs modeled as ellipsoids whose
    %   geometry follows lung_ellipsoid_waveform.
    %
    %   Constructor:
    %       obj = BreathingLung(t_s, pulmonaryOpts, intensity, shapeParameters)
    %
    %   Inputs:
    %       t_s           - time [s]
    %       f_bpm         - breathing frequency [breaths/min]
    %       VT_L          - tidal volume [L]
    %       Vres_L        - residual volume [L]
    %       Vbase_L       - baseline lung volume [L]
    %       bellyFrac     - belly-breathing fraction [0,1]
    %       inspFrac      - inspiratory fraction of cycle (0,1)
    %       pulmonaryOpts.lungSeparation_mm - optional spacing beyond the
    %                      lung radius [mm]
    %       intensity     - shape intensity (handled by parent)
    %       shapeParameters - optional pose struct for the composite
    %
    %   Notes:
    %       - Left and right lungs are created as AnalyticalEllipsoid3D
    %         components and supplied to the CompositeAnalyticalShape3D
    %         constructor as additive components.

    properties (Access = protected)
        rightLung
        leftLung
    end

    methods
        function obj = BreathingLung(t_s, pulmonaryOpts, intensity, shapeParameters)
            if nargin < 3 || isempty(intensity)
                intensity = 1;
            end

            if nargin < 4 || isempty(shapeParameters)
                shapeParameters = BreathingLung.defaultLungParameters();
            else
                shapeParameters = AnalyticalShape3D.ensurePoseFields(shapeParameters);
            end

            waveformHandle = @() BreathingLung.buildLungWaveformParameters(t_s, pulmonaryOpts);
            rightParamsHandle = @() BreathingLung.extractLungParameters(waveformHandle, 'right');
            leftParamsHandle = @() BreathingLung.extractLungParameters(waveformHandle, 'left');

            rightLung = AnalyticalEllipsoid3D([], rightParamsHandle);
            leftLung = AnalyticalEllipsoid3D([], leftParamsHandle);

            obj@CompositeAnalyticalShape3D([leftLung, rightLung], ...
                AnalyticalShape3D.empty, intensity, shapeParameters);

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

    methods (Static, Access = private)
        function [leftParams, rightParams] = buildLungWaveformParameters(t_s, pulmonaryOpts)
            [lungRadius_mm, lungHeight_mm] = lung_ellipsoid_waveform(t_s, pulmonaryOpts);

            if isfield(pulmonaryOpts, 'lungSeparation_mm') && ...
                    ~isempty(pulmonaryOpts.lungSeparation_mm)
                lungSeparation_mm = pulmonaryOpts.lungSeparation_mm;
            else
                lungSeparation_mm = 0;
            end

            if isscalar(lungSeparation_mm)
                lungSeparation_mm = repmat(lungSeparation_mm, size(lungRadius_mm));
            elseif ~isequal(size(lungSeparation_mm), size(lungRadius_mm))
                error('BreathingLung:SeparationSizeMismatch', ...
                    ['pulmonaryOpts.lungSeparation_mm must be scalar or have ' ...
                    'the same size as the lung waveform.']);
            end

            lungPosition_mm = lungRadius_mm + lungSeparation_mm;
            lungParams = struct('a_mm', lungRadius_mm, ...
                'b_mm', lungRadius_mm, ...
                'c_mm', lungHeight_mm, ...
                'R_mm', lungRadius_mm, ...
                'H_mm', lungHeight_mm, ...
                'lungPosition_mm', lungPosition_mm);

            rightParams = lungParams;
            rightParams.pose = struct('center', struct('x_mm', lungPosition_mm(:), ...
                'y_mm', zeros(numel(lungPosition_mm), 1), ...
                'z_mm', zeros(numel(lungPosition_mm), 1)), ...
                'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);

            leftParams = lungParams;
            leftParams.pose = struct('center', struct('x_mm', -lungPosition_mm(:), ...
                'y_mm', zeros(numel(lungPosition_mm), 1), ...
                'z_mm', zeros(numel(lungPosition_mm), 1)), ...
                'roll_deg', 0, 'pitch_deg', 0, 'yaw_deg', 0);
        end

        function params = extractLungParameters(waveformHandle, side)
            [leftParams, rightParams] = waveformHandle();
            switch lower(side)
                case 'left'
                    params = leftParams;
                case 'right'
                    params = rightParams;
                otherwise
                    error('BreathingLung:InvalidSide', 'Side must be ''left'' or ''right''.');
            end
        end

        function params = defaultLungParameters()
            params = AnalyticalShape3D.ensurePoseFields(struct());
        end
    end
end
