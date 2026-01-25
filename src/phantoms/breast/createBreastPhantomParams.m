function params = createBreastPhantomParams()
    % createBreastPhantomParams  Build default parameters for BreastPhantom.
    %
    %   params = createBreastPhantomParams() returns a struct of geometry,
    %   contrast, and timing parameters used by BreastPhantom to construct
    %   the thoracic/breast analytic phantom.

    %% Create heart
    params.cardiacOpts = struct('HR_bpm', 70/10.66, ...
        'EDV_ml', 150, ...
        'ESV_ml', 75, ...
        'systFrac', 0.35, ...
        'q_ED', 50/27, ...
        'GLS_peak', -0.20, ...
        'GCS_peak', -0.25);
    params.heartIntensity = 1;
    params.heartWallThickness_mm = 8;

    %% Lungs
    params.pulmonaryOpts = struct('f_bpm', 12/10.66, ...
        'VT_L', 0.4, ...
        'Vres_L', 0.8, ...
        'Vbase_L', 1.5, ...
        'bellyFrac', 0.6, ...
        'inspFrac', 1/3, ...
        'GCS_peak', 0.4);
    params.lungIntensity = 0.1;

    %% Thorax
    params.tissueGap_lr_mm = 30;
    params.phantomDepth_mm = 300;
    params.fatIntensity = 2;
    params.tissueIntensity = 0.5;

    %% Breast
    params.breast_gap_mm = 60;
    params.breast_radius_mm = 60;
    params.breast_depth_mm = 125;
    params.breastIntensity = 0.5;

    %% Vessel
    params.vesselDiameter_mm = 5;
    params.vesselRadius_mm = 0.5 * params.vesselDiameter_mm;
    params.totalVesselLength_mm = 100;
    params.vesselSegmentCount = 3;
    params.breastRollPitchYaw = [0, 90, 90];
    params.enhancedIntensity = 2.5;
    params.unenhancedIntensity = 0.3;
    params.breastVesselVelocity_cm_s = 1/10.66;
    params.startInjectionTime_s = 30 * 10.66;
    % vesselSegments define connected vessel segments that fill sequentially.
    params.vesselSegments = buildDefaultVesselSegments(params);
end

function segments = buildDefaultVesselSegments(params)
    segmentLength_mm = params.totalVesselLength_mm / params.vesselSegmentCount;
    startZ_mm = -0.5 * params.totalVesselLength_mm + 0.5 * segmentLength_mm;
    segments = repmat(struct('length_mm', segmentLength_mm, ...
        'radius_mm', params.vesselRadius_mm, ...
        'pose', struct()), 1, params.vesselSegmentCount);

    for idx = 1:params.vesselSegmentCount
        centerZ_mm = startZ_mm + (idx - 1) * segmentLength_mm;
        segments(idx).pose = createPoseStruct(0, 0, centerZ_mm, ...
            params.breastRollPitchYaw(1), ...
            params.breastRollPitchYaw(2), ...
            params.breastRollPitchYaw(3));
    end
end

function pose = createPoseStruct(x_mm, y_mm, z_mm, roll, pitch, yaw)
    pose = struct('center', struct('x_mm', x_mm, ...
        'y_mm', y_mm, ...
        'z_mm', z_mm), ...
        'roll_deg', roll, ...
        'pitch_deg', pitch, ...
        'yaw_deg', yaw);
end
