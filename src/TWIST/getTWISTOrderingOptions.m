function orderingOptions = getTWISTOrderingOptions(twistParameters)
% getTWISTOrderingOptions  Normalize TWIST phase-encode ordering options.
%   orderingOptions = getTWISTOrderingOptions(twistParameters) returns a
%   struct with validated fields:
%       radialBinWidthMode : "max" or "min"
%           Selects whether radial shells are quantized using the larger or
%           smaller phase/slice k-space pixel size [cycles / distance].
%       bSubsetAssignment  : "contiguous" or "interleaved"
%           Controls how the ordered region-B list is split into Bj subsets.
%
%   Missing fields fall back to Siemens-like defaults that emphasize
%   shell-wise angular ordering with contiguous Bj continuation across
%   adjacent frames.

    arguments
        twistParameters (1,1) struct
    end

    orderingOptions = struct();

    if isfield(twistParameters, "radialBinWidthMode")
        orderingOptions.radialBinWidthMode = string(twistParameters.radialBinWidthMode);
    else
        orderingOptions.radialBinWidthMode = "max";
    end

    orderingOptions.radialBinWidthMode = lower(orderingOptions.radialBinWidthMode);
    validBinModes = ["max", "min"];
    if ~ismember(orderingOptions.radialBinWidthMode, validBinModes)
        error("TWIST:getTWISTOrderingOptions:InvalidRadialBinWidthMode", ...
            "TWIST.radialBinWidthMode must be one of %s.", ...
            char(strjoin(validBinModes, ", ")));
    end

    if isfield(twistParameters, "bSubsetAssignment")
        orderingOptions.bSubsetAssignment = string(twistParameters.bSubsetAssignment);
    else
        orderingOptions.bSubsetAssignment = "contiguous";
    end

    orderingOptions.bSubsetAssignment = lower(orderingOptions.bSubsetAssignment);
    validAssignments = ["contiguous", "interleaved"];
    if ~ismember(orderingOptions.bSubsetAssignment, validAssignments)
        error("TWIST:getTWISTOrderingOptions:InvalidBSubsetAssignment", ...
            "TWIST.bSubsetAssignment must be one of %s.", ...
            char(strjoin(validAssignments, ", ")));
    end
end
