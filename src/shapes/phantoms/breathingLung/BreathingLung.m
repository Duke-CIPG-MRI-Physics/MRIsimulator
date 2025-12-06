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
        R_mm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        H_mm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end

    methods
        function obj = BreathingLung(provider, lungSeparation_mm, intensity, center, rollPitchYaw)
            arguments
                provider (1,1) PhantomContext
                lungSeparation_mm (1,:) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                intensity double {mustBeFinite} = 0.1
                center (1,3) double {mustBeFinite} = [0, 0, 0]
                rollPitchYaw (1,3) double {mustBeFinite} = [0 0 0]
            end

            R_mm = provider.getLungRadiusMm();
            H_mm = provider.getLungHeightMm();

            lungPosition_mm = R_mm + lungSeparation_mm;

            rightCenter = [lungPosition_mm(:), zeros(size(lungPosition_mm(:))), zeros(size(lungPosition_mm(:)))];
            leftCenter = [-lungPosition_mm(:), zeros(size(lungPosition_mm(:))), zeros(size(lungPosition_mm(:)))];


            rightLung = AnalyticalEllipsoid3D(R_mm, R_mm, H_mm, [], rightCenter, [0, 0, 0]);
            leftLung = AnalyticalEllipsoid3D(R_mm, R_mm, H_mm, [], leftCenter, [0, 0, 0]);

            rightLung.setTimeSamples(t_s);
            leftLung.setTimeSamples(t_s);

            obj@SharedIntensityShapeGroup3D([leftLung, rightLung], ...
                AnalyticalShape3D.empty, intensity, center, rollPitchYaw);

            obj.R_mm = R_mm;
            obj.H_mm = H_mm;
        end

        function lungRadius = getLungRadiusMm(obj)
            lungRadius = obj.R_mm;
        end

        function lungHeight = getLungHeightMm(obj)
            lungHeight = obj.H_mm;
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
end
