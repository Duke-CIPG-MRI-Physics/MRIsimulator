classdef ShapeProviderTests < matlab.unittest.TestCase
    % ShapeProviderTests
    %   Verify that dependent geometry relationships are preserved when
    %   shapes query their providers directly at evaluation time.

    methods (TestClassSetup)
        function addSourcePaths(~)
            baseDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(baseDir));
        end
    end

    methods (Test)
        function lungsFollowHeartExpansion(testCase)
            t_s = linspace(0, 1, 6);
            heart = BeatingHeart(t_s, 70 * ones(size(t_s)), ...
                140 * ones(size(t_s)), 70 * ones(size(t_s)), 1);

            lungSpacing_mm = heart.getA() + 8;  % mimic BreastPhantom coupling
            lungs = BreathingLung(t_s, 12 * ones(size(t_s)), ...
                0.4 * ones(size(t_s)), 0.8 * ones(size(t_s)), ...
                1.5 * ones(size(t_s)), 0.6 * ones(size(t_s)), ...
                0.4 * ones(size(t_s)), lungSpacing_mm, 0.1);

            centers = lungs.getLungCenters();
            expectedOffset = lungs.getLungRadiusMm().' + lungSpacing_mm(:);

            testCase.verifyEqual(centers.right(:,1), expectedOffset, 'AbsTol', 1e-12);
            testCase.verifyEqual(centers.left(:,1), -expectedOffset, 'AbsTol', 1e-12);
        end

        function vesselContrastTracksVolume(testCase)
            t_s = (0:4).';
            radius_mm = 2.5;
            totalLength_mm = 50;
            maxVolume_mm3 = pi * radius_mm^2 * totalLength_mm;
            contrastVolume_mm3 = linspace(0, maxVolume_mm3, numel(t_s)).';

            vessel = EnhancingVessel(t_s, totalLength_mm, 2.5, 0.4, ...
                radius_mm, contrastVolume_mm3);

            [enhanced, unenhanced] = vessel.getVessels();
            expectedEnhanced = computeContrastWashIn(t_s, radius_mm, contrastVolume_mm3);

            testCase.verifyEqual(enhanced.getLength(), expectedEnhanced, 'AbsTol', 1e-12);
            testCase.verifyEqual(unenhanced.getLength(), ...
                totalLength_mm - expectedEnhanced, 'AbsTol', 1e-12);
        end
    end
end
