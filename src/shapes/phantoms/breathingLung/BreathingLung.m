classdef BreathingLung < SharedIntensityShapeGroup3D
    % BreathingLung
    %   Composite of two breathing lungs modeled as ellipsoids whose
    %   geometry follows computeBreathingMotionEllipsoid.
    %
    %   Constructor:
    %       obj = BreathingLung(provider, lungSeparation_mm, intensity, ...)
    %
    %   Inputs:
    %       provider         - PhantomContext supplying lung geometry
    %       lungSeparation_mm- spacing between lungs [mm]
    %       intensity        - shape intensity (handled by parent)
    %       center           - composite center [mm]
    %       rollPitchYaw     - optional composite orientation [deg]

    properties (Access = private)
        radiusFcn (1,1) function_handle
        heightFcn (1,1) function_handle
        separationFcn (1,1) function_handle
        time_s (1,:) double {mustBeReal, mustBeFinite}
    end

    methods
        function obj = BreathingLung(provider, lungSeparation_mm, intensity, center, rollPitchYaw)
            arguments
                provider (1,1) PhantomContext
                lungSeparation_mm
                intensity double {mustBeFinite} = 0.1
                center (1,3) double {mustBeFinite} = [0, 0, 0]
                rollPitchYaw (1,3) double {mustBeFinite} = [0 0 0]
            end

            obj.radiusFcn = provider.getLungRadiusMm();
            obj.heightFcn = provider.getLungHeightMm();
            obj.separationFcn = obj.normalizeSeparation(lungSeparation_mm);
            obj.time_s = provider.getTime();

            lungPosition_mm = obj.radiusFcn(obj.time_s) + obj.separationFcn(obj.time_s);

            rightCenter = [lungPosition_mm(:), zeros(size(lungPosition_mm(:))), zeros(size(lungPosition_mm(:)))];
            leftCenter = [-lungPosition_mm(:), zeros(size(lungPosition_mm(:))), zeros(size(lungPosition_mm(:)))];

            rightLung = AnalyticalEllipsoid3D([], [], [], [], rightCenter, [0, 0, 0]);
            leftLung = AnalyticalEllipsoid3D([], [], [], [], leftCenter, [0, 0, 0]);

            opts = struct('Cache', false);
            rightLung.setA(obj.radiusFcn, opts);
            rightLung.setB(obj.radiusFcn, opts);
            rightLung.setC(obj.heightFcn, opts);

            leftLung.setA(obj.radiusFcn, opts);
            leftLung.setB(obj.radiusFcn, opts);
            leftLung.setC(obj.heightFcn, opts);

            rightLung.setTimeSamples(obj.time_s);
            leftLung.setTimeSamples(obj.time_s);

            obj@SharedIntensityShapeGroup3D([leftLung, rightLung], ...
                AnalyticalShape3D.empty, intensity, center, rollPitchYaw);
        end

        function lungRadius = getLungRadiusMm(obj, t_s)
            lungRadius = obj.radiusFcn(obj.normalizeTime(t_s));
        end

        function lungHeight = getLungHeightMm(obj, t_s)
            lungHeight = obj.heightFcn(obj.normalizeTime(t_s));
        end

        function centers = getLungCenters(obj)
            % getLungCenters
            %   Return WORLD centers for the left and right lung components.

            if numel(obj.additiveComponents) < 2
                error('BreathingLung:MissingComponents', ...
                    'Expected left and right lung components to be present.');
            end

            centers.left = obj.additiveComponents(1).getCenter();
            centers.right = obj.additiveComponents(2).getCenter();
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

        function separationFcn = normalizeSeparation(~, lungSeparation_mm)
            if isa(lungSeparation_mm, 'function_handle')
                separationFcn = @(t) lungSeparation_mm(t);
            else
                validateattributes(lungSeparation_mm, {'double'}, {'nonnegative', 'finite'});
                separationFcn = @(t) lungSeparation_mm .* ones(size(t));
            end
        end
    end
end
