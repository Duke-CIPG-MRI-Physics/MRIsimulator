function segments = createBreastVesselSegments(segmentLengths_mm, radius_mm, vesselOffsetRl_mm, ...
    chainCenter_mm, chainRollPitchYaw_deg)
    % createBreastVesselSegments  Build contiguous vessel segments for BreastPhantom.
    %
    %   segments = createBreastVesselSegments(segmentLengths_mm, radius_mm, vesselOffsetRl_mm)
    %   returns a struct array of connected vessel segments offset along the
    %   right-left (x) axis with no rotation.
    %
    %   segments = createBreastVesselSegments(segmentLengths_mm, radius_mm, ...
    %       vesselOffsetRl_mm, chainCenter_mm, chainRollPitchYaw_deg) sets the chain center
    %   [x y z] in mm and roll/pitch/yaw in degrees for the vessel axis. Even
    %   segments are rotated 180 degrees about the y-axis (pitch) so that
    %   every other vessel flips orientation.
    %
    %   Inputs:
    %       segmentLengths_mm      - Row vector of segment lengths [mm].
    %       radius_mm              - Scalar vessel radius [mm].
    %       vesselOffsetRl_mm      - Scalar right-left separation [mm].
    %       chainCenter_mm         - Optional [x y z] chain center [mm].
    %       chainRollPitchYaw_deg  - Optional [roll pitch yaw] in degrees.
    %
    %   Output:
    %       segments - Struct array with fields length_mm, radius_mm, and pose.
    %
    %   Example:
    %       segments = createBreastVesselSegments([30 40 30], 2.5, 20, [0 0 0], [0 90 90]);

    if nargin < 3
        vesselOffsetRl_mm = 0;
    end
    if nargin < 4
        chainCenter_mm = [0, 0, 0];
    end
    if nargin < 5
        chainRollPitchYaw_deg = [0, 0, 0];
    end

    validateattributes(segmentLengths_mm, {'numeric'}, {'row', 'real', 'finite', 'positive'}, ...
        mfilename, 'segmentLengths_mm');
    validateattributes(radius_mm, {'numeric'}, {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'radius_mm');
    validateattributes(vesselOffsetRl_mm, {'numeric'}, {'real', 'finite', 'scalar'}, ...
        mfilename, 'vesselOffsetRl_mm');
    validateattributes(chainCenter_mm, {'numeric'}, {'real', 'finite', 'numel', 3}, ...
        mfilename, 'chainCenter_mm');
    validateattributes(chainRollPitchYaw_deg, {'numeric'}, {'real', 'finite', 'numel', 3}, ...
        mfilename, 'chainRollPitchYaw_deg');

    segmentCount = numel(segmentLengths_mm);
    centerOffsetsRl_mm = ((1:segmentCount) - (segmentCount + 1) / 2) * vesselOffsetRl_mm;

    segments = repmat(struct('length_mm', 0, 'radius_mm', radius_mm, 'pose', struct()), ...
        1, numel(segmentLengths_mm));

    for idx = 1:numel(segmentLengths_mm)
        center_mm = chainCenter_mm(:)';
        center_mm(1) = center_mm(1) + centerOffsetsRl_mm(idx);
        segmentPitch_deg = chainRollPitchYaw_deg(2);
        if mod(idx, 2) == 0
            segmentPitch_deg = segmentPitch_deg + 180;
        end
        segments(idx).length_mm = segmentLengths_mm(idx);
        segments(idx).radius_mm = radius_mm;
        segments(idx).pose = struct('center', struct('x_mm', center_mm(1), ...
            'y_mm', center_mm(2), ...
            'z_mm', center_mm(3)), ...
            'roll_deg', chainRollPitchYaw_deg(1), ...
            'pitch_deg', segmentPitch_deg, ...
            'yaw_deg', chainRollPitchYaw_deg(3));
    end
end
