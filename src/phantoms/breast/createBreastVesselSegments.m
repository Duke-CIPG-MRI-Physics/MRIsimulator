function segments = createBreastVesselSegments(segmentLengths_mm, radius_mm, chainCenter_mm, chainRollPitchYaw_deg)
    % createBreastVesselSegments  Build contiguous vessel segments for BreastPhantom.
    %
    %   segments = createBreastVesselSegments(segmentLengths_mm, radius_mm)
    %   returns a struct array of connected vessel segments centered about the
    %   origin with no rotation.
    %
    %   segments = createBreastVesselSegments(segmentLengths_mm, radius_mm, ...
    %       chainCenter_mm, chainRollPitchYaw_deg) sets the chain center
    %   [x y z] in mm and roll/pitch/yaw in degrees for the vessel axis.
    %
    %   Inputs:
    %       segmentLengths_mm      - Row vector of segment lengths [mm].
    %       radius_mm              - Scalar vessel radius [mm].
    %       chainCenter_mm         - Optional [x y z] chain center [mm].
    %       chainRollPitchYaw_deg  - Optional [roll pitch yaw] in degrees.
    %
    %   Output:
    %       segments - Struct array with fields length_mm, radius_mm, and pose.
    %
    %   Example:
    %       segments = createBreastVesselSegments([30 40 30], 2.5, [0 0 0], [0 90 90]);

    if nargin < 3
        chainCenter_mm = [0, 0, 0];
    end
    if nargin < 4
        chainRollPitchYaw_deg = [0, 0, 0];
    end

    validateattributes(segmentLengths_mm, {'numeric'}, {'row', 'real', 'finite', 'positive'}, ...
        mfilename, 'segmentLengths_mm');
    validateattributes(radius_mm, {'numeric'}, {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'radius_mm');
    validateattributes(chainCenter_mm, {'numeric'}, {'real', 'finite', 'numel', 3}, ...
        mfilename, 'chainCenter_mm');
    validateattributes(chainRollPitchYaw_deg, {'numeric'}, {'real', 'finite', 'numel', 3}, ...
        mfilename, 'chainRollPitchYaw_deg');

    axisUnit = axisUnitFromRollPitchYaw(chainRollPitchYaw_deg);
    totalLength_mm = sum(segmentLengths_mm);
    startOffset_mm = -0.5 * totalLength_mm;
    cumulativeLengths_mm = [0, cumsum(segmentLengths_mm(1:end-1))];
    centerOffsets_mm = startOffset_mm + cumulativeLengths_mm + 0.5 .* segmentLengths_mm;

    segments = repmat(struct('length_mm', 0, 'radius_mm', radius_mm, 'pose', struct()), ...
        1, numel(segmentLengths_mm));

    for idx = 1:numel(segmentLengths_mm)
        center_mm = chainCenter_mm(:)' + axisUnit(:)' .* centerOffsets_mm(idx);
        segments(idx).length_mm = segmentLengths_mm(idx);
        segments(idx).radius_mm = radius_mm;
        segments(idx).pose = struct('center', struct('x_mm', center_mm(1), ...
            'y_mm', center_mm(2), ...
            'z_mm', center_mm(3)), ...
            'roll_deg', chainRollPitchYaw_deg(1), ...
            'pitch_deg', chainRollPitchYaw_deg(2), ...
            'yaw_deg', chainRollPitchYaw_deg(3));
    end
end

function axisUnit = axisUnitFromRollPitchYaw(chainRollPitchYaw_deg)
    roll_rad = deg2rad(chainRollPitchYaw_deg(1));
    pitch_rad = deg2rad(chainRollPitchYaw_deg(2));
    yaw_rad = deg2rad(chainRollPitchYaw_deg(3));

    cr = cos(roll_rad);
    sr = sin(roll_rad);
    cp = cos(pitch_rad);
    sp = sin(pitch_rad);
    cy = cos(yaw_rad);
    sy = sin(yaw_rad);

    axisUnit = [cy .* sp .* cr + sy .* sr, ...
        sy .* sp .* cr - cy .* sr, ...
        cp .* cr];
end
