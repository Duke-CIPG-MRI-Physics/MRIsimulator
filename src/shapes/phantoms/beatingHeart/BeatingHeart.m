classdef BeatingHeart < AnalyticalEllipsoid3D
    % BeatingHeart
    %   Time-varying ellipsoidal heart model driven by cardiac_ellipsoid_waveform.
    %
    %   Constructor:
    %       heart = BeatingHeart(provider, intensity, center, rollPitchYaw)
    %
    %   Inputs:
    %       provider      : PhantomContext supplying precomputed radii
    %       intensity     : (optional) shape intensity scaling
    %       center        : (optional) WORLD center [mm]
    %       rollPitchYaw  : (optional) [roll pitch yaw] in deg

    methods
        function obj = BeatingHeart(provider, intensity, center, rollPitchYaw)
            arguments
                provider (1,1) PhantomContext
                intensity double {mustBeReal, mustBeFinite} = 1
                center double {mustBeReal, mustBeFinite} = [0 0 0]
                rollPitchYaw (1,3) double {mustBeReal, mustBeFinite} = [0 0 0]
            end

            [aFcn, bFcn, cFcn] = provider.getHeartRadiiMm();

            obj@AnalyticalEllipsoid3D([], [], [], intensity, center, ...
                rollPitchYaw);

            opts = struct('Cache', false);
            obj.setA(aFcn, opts);
            obj.setB(bFcn, opts);
            obj.setC(cFcn, opts);
            obj.setTimeSamples(provider.getTime());
        end
    end
end
