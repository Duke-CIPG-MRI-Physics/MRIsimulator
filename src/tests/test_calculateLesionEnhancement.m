function test_calculateLesionEnhancement()
% test_calculateLesionEnhancement  Deterministic checks for lesion kinetics.
%
%   test_calculateLesionEnhancement() validates that empirical curve presets
%   produce expected late-phase ranking for persistent, plateau, and washout
%   lesion behavior.

    params = createBreastPhantomParams();
    params.startInjectionTime_s = 0;
    params.lesionArrivalDelay_s = 8;
    params.lesionPeakEnhancement = 1.6;
    params.lesionBaselineDeltaIntensity = 0;

    t_s = linspace(0, 300, 400);

    persistentCurve = calculateLesionEnhancement(t_s, params, "fast", "persistent");
    plateauCurve = calculateLesionEnhancement(t_s, params, "fast", "plateau");
    washoutCurve = calculateLesionEnhancement(t_s, params, "fast", "washout");

    lateMask = t_s > 220;
    persistentLate = mean(persistentCurve(lateMask));
    plateauLate = mean(plateauCurve(lateMask));
    washoutLate = mean(washoutCurve(lateMask));

    assert(persistentLate > plateauLate, ...
        'Expected persistent late enhancement to exceed plateau behavior.');
    assert(plateauLate > washoutLate, ...
        'Expected plateau late enhancement to exceed washout behavior.');

    assert(all(persistentCurve >= 0) && all(plateauCurve >= 0) && all(washoutCurve >= 0), ...
        'Expected non-negative lesion intensity outputs.');
end
