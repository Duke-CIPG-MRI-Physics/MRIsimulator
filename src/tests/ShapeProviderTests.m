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
            t_s = linspace(0, 1, 6).';

            heartA = @(t) 70 * ones(size(t));
            heartB = @(t) 140 * ones(size(t));
            heartC = @(t) 70 * ones(size(t));

            heart = BeatingHeart(heartA, heartB, heartC, t_s.');

            lungSpacing_mm = @(t) heartA(t) + 8;  % mimic BreastPhantom coupling
            radiusFcn = @(t) 12 * ones(size(t));
            heightFcn = @(t) 0.4 * ones(size(t));
            lungs = BreathingLung(radiusFcn, heightFcn, lungSpacing_mm, t_s.', 0.1);

            centers = lungs.getLungCenters();
            expectedOffset = lungs.getLungRadiusMm(t_s.').' + lungSpacing_mm(t_s.').';

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
