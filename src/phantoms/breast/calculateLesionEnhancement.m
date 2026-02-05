function lesionIntensity = calculateLesionEnhancement(t_s, breastPhantomParams, washinType, washoutType, kineticOverrides)
% calculateLesionEnhancement  Compute simplified lesion enhancement kinetics.
%
%   lesionIntensity = calculateLesionEnhancement(t_s, breastPhantomParams, ...
%       washinType, washoutType) returns lesion-only enhancement intensity
%   [a.u.] over time t_s [s].
%
%   lesionIntensity = calculateLesionEnhancement(..., kineticOverrides)
%   allows optional override of kinetic parameters:
%       kWashin_per_s      : wash-in rate [1/s]
%       kWashout_per_s     : washout rate [1/s]
%       timeToPeak_s       : time from arrival to transition [s]
%       latePhaseFraction  : late level as fraction of peak [0..1.5]
%
%   washinType controls wash-in speed: "slow", "medium", or "fast".
%   washoutType controls late-phase shape:
%       "persistent" (type I), "plateau" (type II), "washout" (type III).
%
%   Example
%   -------
%   params = createBreastPhantomParams();
%   t_s = linspace(0, 300, 200);
%   lesionIntensity = calculateLesionEnhancement(t_s, params, "fast", "washout");

    arguments
        t_s {mustBeNumeric, mustBeReal, mustBeFinite}
        breastPhantomParams struct
        washinType {mustBeTextScalar} = "medium"
        washoutType {mustBeTextScalar} = "plateau"
        kineticOverrides struct = struct()
    end

    %% Input validation
    validateattributes(t_s, {'numeric'}, {'nonempty'}, mfilename, 't_s');
    washinType = lower(string(washinType));
    washoutType = lower(string(washoutType));

    %% Resolve model parameters
    modelParams = getDefaultKineticParams(washinType, washoutType);
    modelParams = mergeStructs(modelParams, kineticOverrides);

    startInjectionTime_s = getFieldWithDefault(breastPhantomParams, 'startInjectionTime_s', 0);
    lesionArrivalDelay_s = getFieldWithDefault(breastPhantomParams, 'lesionArrivalDelay_s', 8);
    lesionPeakEnhancement = getFieldWithDefault(breastPhantomParams, 'lesionPeakEnhancement', 1.6);
    lesionBaselineDeltaIntensity = getFieldWithDefault(breastPhantomParams, 'lesionBaselineDeltaIntensity', 0);

    validateattributes(lesionPeakEnhancement, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}, ...
        mfilename, 'lesionPeakEnhancement');

    %% Core logic
    tArrival_s = startInjectionTime_s + lesionArrivalDelay_s;
    tau_s = max(t_s - tArrival_s, 0);

    peakDelta = max(lesionPeakEnhancement - lesionBaselineDeltaIntensity, 0);
    washinFraction = 1 - exp(-modelParams.kWashin_per_s .* tau_s);
    candidateSignal = lesionBaselineDeltaIntensity + peakDelta .* washinFraction;

    tauLate_s = max(tau_s - modelParams.timeToPeak_s, 0);
    peakSignal = lesionBaselineDeltaIntensity + peakDelta .* ...
        (1 - exp(-modelParams.kWashin_per_s .* modelParams.timeToPeak_s));
    lateTargetSignal = lesionBaselineDeltaIntensity + modelParams.latePhaseFraction .* peakDelta;

    washoutWeight = 1 - exp(-modelParams.kWashout_per_s .* tauLate_s);
    lateSignal = peakSignal - (peakSignal - lateTargetSignal) .* washoutWeight;

    lesionIntensity = candidateSignal;
    isLatePhase = tau_s > modelParams.timeToPeak_s;
    lesionIntensity(isLatePhase) = lateSignal(isLatePhase);

    %% Output assembly
    lesionIntensity = max(lesionIntensity, 0);
end

function modelParams = getDefaultKineticParams(washinType, washoutType)
    switch washinType
        case "slow"
            kWashin_per_s = 0.020;
            timeToPeak_s = 140;
        case "medium"
            kWashin_per_s = 0.040;
            timeToPeak_s = 90;
        case "fast"
            kWashin_per_s = 0.070;
            timeToPeak_s = 55;
        otherwise
            error('calculateLesionEnhancement:InvalidWashinType', ...
                'Unknown washinType "%s". Use slow, medium, or fast.', washinType);
    end

    switch washoutType
        case "persistent"
            latePhaseFraction = 1.10;
            kWashout_per_s = 0.004;
        case "plateau"
            latePhaseFraction = 0.95;
            kWashout_per_s = 0.008;
        case "washout"
            latePhaseFraction = 0.60;
            kWashout_per_s = 0.018;
        otherwise
            error('calculateLesionEnhancement:InvalidWashoutType', ...
                'Unknown washoutType "%s". Use persistent, plateau, or washout.', washoutType);
    end

    modelParams = struct( ...
        'kWashin_per_s', kWashin_per_s, ...
        'kWashout_per_s', kWashout_per_s, ...
        'timeToPeak_s', timeToPeak_s, ...
        'latePhaseFraction', latePhaseFraction);
end

function value = getFieldWithDefault(inputStruct, fieldName, defaultValue)
    if isfield(inputStruct, fieldName) && ~isempty(inputStruct.(fieldName))
        value = inputStruct.(fieldName);
    else
        value = defaultValue;
    end
end

function outputStruct = mergeStructs(defaultStruct, overrideStruct)
    outputStruct = defaultStruct;
    overrideFields = fieldnames(overrideStruct);
    for iField = 1:numel(overrideFields)
        fieldName = overrideFields{iField};
        outputStruct.(fieldName) = overrideStruct.(fieldName);
    end
end
