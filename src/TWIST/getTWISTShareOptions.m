function shareOptions = getTWISTShareOptions(twistParameters)
% getTWISTShareOptions  Normalize TWIST view-sharing options.
%   shareOptions = getTWISTShareOptions(twistParameters) returns a struct
%   with validated fields:
%       mode       : "forward", "reverse", or "symmetric"
%       method     : "single_anchor" or "dual_anchor"
%       tieBreaker : "future" or "past" for symmetric equidistant cases
%
%   Missing fields fall back to Siemens-like defaults:
%       forward + single_anchor for legacy behavior
%       symmetric + dual_anchor when symmetric mode is requested

arguments
    twistParameters (1,1) struct
end

shareOptions = struct;

if isfield(twistParameters, "shareMode")
    shareOptions.mode = string(twistParameters.shareMode);
else
    shareOptions.mode = "forward";
end

shareOptions.mode = lower(shareOptions.mode);
validModes = ["forward", "reverse", "symmetric"];
if ~ismember(shareOptions.mode, validModes)
    error("TWIST:getTWISTShareOptions:InvalidShareMode", ...
        "TWIST.shareMode must be one of %s.", char(strjoin(validModes, ", ")));
end

if isfield(twistParameters, "shareMethod")
    shareOptions.method = string(twistParameters.shareMethod);
elseif shareOptions.mode == "symmetric"
    shareOptions.method = "dual_anchor";
else
    shareOptions.method = "single_anchor";
end

shareOptions.method = lower(shareOptions.method);
validMethods = ["single_anchor", "dual_anchor"];
if ~ismember(shareOptions.method, validMethods)
    error("TWIST:getTWISTShareOptions:InvalidShareMethod", ...
        "TWIST.shareMethod must be one of %s.", char(strjoin(validMethods, ", ")));
end

if shareOptions.mode == "symmetric" && shareOptions.method ~= "dual_anchor"
    error("TWIST:getTWISTShareOptions:UnsupportedSymmetricMethod", ...
        "Symmetric TWIST sharing currently requires TWIST.shareMethod = ""dual_anchor"".");
end

if isfield(twistParameters, "shareTieBreaker")
    shareOptions.tieBreaker = string(twistParameters.shareTieBreaker);
else
    % Siemens flexible view sharing descriptions prefer backward sharing
    % when two neighboring measurements are equally close in time.
    shareOptions.tieBreaker = "future";
end

shareOptions.tieBreaker = lower(shareOptions.tieBreaker);
validTieBreakers = ["future", "past"];
if ~ismember(shareOptions.tieBreaker, validTieBreakers)
    error("TWIST:getTWISTShareOptions:InvalidTieBreaker", ...
        "TWIST.shareTieBreaker must be one of %s.", char(strjoin(validTieBreakers, ", ")));
end
end
