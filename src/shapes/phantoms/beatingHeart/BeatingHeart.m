classdef BeatingHeart < AnalyticalEllipsoid3D
    % BeatingHeart
    %   Time-varying ellipsoidal heart model driven by cardiac_ellipsoid_waveform.
    %
    %   Constructor:
    %       heart = BeatingHeart(aFcn, bFcn, cFcn, time_s, intensity, center, rollPitchYaw)
    %
    %   Inputs:
    %       aFcn, bFcn, cFcn : function handles returning semi-axis lengths [mm]
    %       time_s           : row vector of time samples [s]
    %       intensity        : (optional) shape intensity scaling
    %       center           : (optional) WORLD center [mm]
    %       rollPitchYaw     : (optional) [roll pitch yaw] in deg

    methods
        function obj = BeatingHeart(aFcn, bFcn, cFcn, time_s, intensity, center, rollPitchYaw)
            arguments
                aFcn (1,1) function_handle
                bFcn (1,1) function_handle
                cFcn (1,1) function_handle
                time_s (1,:) double {mustBeReal, mustBeFinite}
                intensity double {mustBeReal, mustBeFinite} = 1
                center double {mustBeReal, mustBeFinite} = [0 0 0]
                rollPitchYaw (1,3) double {mustBeReal, mustBeFinite} = [0 0 0]
            end

            timeRow = double(time_s(:).');

            obj@AnalyticalEllipsoid3D([], [], [], intensity, center, ...
                rollPitchYaw);

            obj.setA(aFcn);
            obj.setB(bFcn);
            obj.setC(cFcn);
            obj.setTimeSamples(timeRow);
        end
    end
end
