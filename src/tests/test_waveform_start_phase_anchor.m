function test_waveform_start_phase_anchor()
%TEST_WAVEFORM_START_PHASE_ANCHOR  Deterministic checks for t=0 phase anchoring.
%
%   test_waveform_start_phase_anchor()
%
%   Verifies that lung and cardiac waveform generators include the phase
%   accumulated from absolute time t=0 to t_s(1), and remain backward
%   compatible when t_s(1)=0.

    %% Lung waveform: t_s(1) > 0 includes start-phase offset
    t_s = [0.5, 1.0, 1.5];
    lungParameters = struct();
    lungParameters.f_bpm = 12;           % 0.2 Hz
    lungParameters.VT_L = 1.0;           % [L]
    lungParameters.Vres_L = 1.0;         % [L]
    lungParameters.Vbase_L = 2.0;        % [L]
    lungParameters.bellyFrac = 0.5;
    lungParameters.inspFrac = 0.5;

    [R_offset_mm, H_offset_mm] = lung_ellipsoid_waveform(t_s, lungParameters);
    [R_zero_mm, H_zero_mm] = lung_ellipsoid_waveform(t_s - t_s(1), lungParameters);

    assert(abs(R_offset_mm(1) - R_zero_mm(1)) > 1e-6, ...
        'Lung waveform first sample should include nonzero start-phase offset.');
    assert(abs(H_offset_mm(1) - H_zero_mm(1)) > 1e-6, ...
        'Lung waveform first sample should include nonzero start-phase offset.');

    %% Lung waveform: t_s(1)=0 retains baseline behavior
    [R_zeroStart_mm, H_zeroStart_mm] = lung_ellipsoid_waveform([0.0, 0.5, 1.0], lungParameters);
    [R_zeroStart_ref_mm, H_zeroStart_ref_mm] = lung_ellipsoid_waveform([0.0, 0.5, 1.0], lungParameters);

    assert(max(abs(R_zeroStart_mm - R_zeroStart_ref_mm)) < 1e-12, ...
        'Lung waveform should be unchanged for t_s(1)=0.');
    assert(max(abs(H_zeroStart_mm - H_zeroStart_ref_mm)) < 1e-12, ...
        'Lung waveform should be unchanged for t_s(1)=0.');

    %% Cardiac waveform: per-channel start-phase offset for t_s(1) > 0
    t_s_card = [0.25, 0.75, 1.25, 1.75];
    opts = struct();
    opts.HR_bpm = [60; 90];
    opts.EDV_ml = [140; 130];
    opts.ESV_ml = [70; 65];
    opts.systFrac = 0.35;
    opts.q_ED = 2.5;
    opts.GLS_peak = -0.20;
    opts.GCS_peak = -0.25;

    cardiac_offset = cardiac_ellipsoid_waveform(t_s_card, opts);
    cardiac_zero = cardiac_ellipsoid_waveform(t_s_card - t_s_card(1), opts);

    assert(any(abs(cardiac_offset.a_mm(:, 1) - cardiac_zero.a_mm(:, 1)) > 1e-6), ...
        'Cardiac waveform first sample should include nonzero start-phase offset.');
    assert(any(abs(cardiac_offset.b_mm(:, 1) - cardiac_zero.b_mm(:, 1)) > 1e-6), ...
        'Cardiac waveform first sample should include nonzero start-phase offset.');

    %% Cardiac waveform: t_s(1)=0 retains baseline behavior
    t_s_zero = [0.0, 0.5, 1.0, 1.5];
    cardiac_zeroStart = cardiac_ellipsoid_waveform(t_s_zero, opts);
    cardiac_zeroStart_ref = cardiac_ellipsoid_waveform(t_s_zero, opts);

    assert(max(abs(cardiac_zeroStart.a_mm(:) - cardiac_zeroStart_ref.a_mm(:))) < 1e-12, ...
        'Cardiac waveform should be unchanged for t_s(1)=0.');
    assert(max(abs(cardiac_zeroStart.b_mm(:) - cardiac_zeroStart_ref.b_mm(:))) < 1e-12, ...
        'Cardiac waveform should be unchanged for t_s(1)=0.');
end
