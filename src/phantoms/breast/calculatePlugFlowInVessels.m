function plugFlowLength_mm = calculatePlugFlowInVessels(t_s, params)
    % calculatePlugFlowInVessels  Compute plug-flow distance through vessels.
    %
    %   plugFlowLength_mm = calculatePlugFlowInVessels(t_s, params) returns
    %   the contrast plug-flow length in mm for each time in t_s, based on
    %   the injection start time and vessel velocity.
    %
    %   Inputs:
    %       t_s    - Time samples [s].
    %       params - Struct with fields:
    %               startInjectionTime_s       [s]
    %               breastVesselVelocity_cm_s  [cm/s]
    %
    %   Output:
    %       plugFlowLength_mm - Contrast plug-flow length [mm].
    %
    %   Example:
    %       params = createBreastPhantomParams();
    %       t_s = linspace(0, 30, 300);
    %       plugFlowLength_mm = calculatePlugFlowInVessels(t_s, params);

    arguments
        t_s {mustBeNumeric, mustBeReal, mustBeFinite}
        params struct
    end

    if ~isfield(params, 'startInjectionTime_s')
        error('calculatePlugFlowInVessels:MissingStartInjectionTime', ...
            'params.startInjectionTime_s is required (seconds).');
    end
    if ~isfield(params, 'breastVesselVelocity_cm_s')
        error('calculatePlugFlowInVessels:MissingVesselVelocity', ...
            'params.breastVesselVelocity_cm_s is required (cm/s).');
    end

    validateattributes(params.startInjectionTime_s, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'params.startInjectionTime_s');
    validateattributes(params.breastVesselVelocity_cm_s, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'params.breastVesselVelocity_cm_s');

    timePostInj_s = max(t_s - params.startInjectionTime_s, 0);
    plugFlowLength_mm = params.breastVesselVelocity_cm_s * 10 .* timePostInj_s;
end
