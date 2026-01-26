function params = createBreastPhantomParams()
    % createBreastPhantomParams  Build default parameters for BreastPhantom.
    %
    %   params = createBreastPhantomParams() returns a struct of geometry,
    %   contrast, and timing parameters used by BreastPhantom to construct
    %   the thoracic/breast analytic phantom.
    %
    %   Vessel segments:
    %       The vessel is modeled as a chain of connected segments that fill
    %       sequentially. Use createBreastVesselSegments to generate a
    %       non-overlapping chain aligned to a center point and roll/pitch/yaw.
    %       By default, segments are offset in the right-left (x) direction
    %       by vessel_offset_rl_mm, with every other segment rotated 180 deg
    %       about the y-axis (pitch). The left chain is additionally rotated
    %       180 deg about the z-axis (yaw) to mirror the right chain.
    %
    %   Example:
    %       params = createBreastPhantomParams();
    %       segmentLengths_mm = [30, 40, 30];
    %       params.vesselSegmentsRight = createBreastVesselSegments( ...
    %           segmentLengths_mm, params.vesselRadius_mm, params.vessel_offset_rl_mm, ...
    %           params.vesselChainCenterRight_mm, params.vesselRollPitchYawRight_deg);
    %       params.vesselSegmentsLeft = createBreastVesselSegments( ...
    %           segmentLengths_mm, params.vesselRadius_mm, params.vessel_offset_rl_mm, ...
    %           params.vesselChainCenterLeft_mm, params.vesselRollPitchYawLeft_deg);

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
    params.vesselSegmentCount = 3;
    params.individualVesselLength_mm = 100 / params.vesselSegmentCount;
    params.totalVesselLength_mm = params.individualVesselLength_mm * params.vesselSegmentCount;
    params.vessel_offset_rl_mm = 20;
    params.breastRollPitchYaw = [0, 90, 90];
    params.enhancedIntensity = 2.5;
    params.unenhancedIntensity = 0.3;
    params.breastVesselVelocity_cm_s = 1/10.66;
    params.startInjectionTime_s = 30 * 10.66;
    params.vesselRollPitchYawRight_deg = params.breastRollPitchYaw;
    params.vesselRollPitchYawLeft_deg = params.breastRollPitchYaw + [0 0 180];
    params.vesselChainCenterRight_mm = [0 0 0];
    params.vesselChainCenterLeft_mm = [0 0 0];
    params.vesselPhantomRollPitchYawRight_deg = [0 0 0];
    params.vesselPhantomRollPitchYawLeft_deg = [0 0 90];
    segmentLengths_mm = repmat(params.individualVesselLength_mm, 1, params.vesselSegmentCount);
    params.vesselSegmentsRight = createBreastVesselSegments( ...
        segmentLengths_mm, params.vesselRadius_mm, params.vessel_offset_rl_mm, ...
        params.vesselChainCenterRight_mm, params.vesselRollPitchYawRight_deg);
    params.vesselSegmentsLeft = createBreastVesselSegments( ...
        segmentLengths_mm, params.vesselRadius_mm, params.vessel_offset_rl_mm, ...
        params.vesselChainCenterLeft_mm, params.vesselRollPitchYawLeft_deg);
    params.vesselSegments = params.vesselSegmentsRight;
end
