clear; clc; close all;
% Dedicated TWIST ordering diagnostic demo.

% --- Inputs ---
pA = 0.05; 
pB = 0.1; 
Matrix_Size_Acquired = [1,100,100];
FOV_acquired = [1,200,100];
R = [1,1];
PF_Factor = [1,1]; 
radialBinWidthMode = "max";

% --- Run TWIST ---
[TWIST_sampling_order] = TWIST( ...
    pA, pB, Matrix_Size_Acquired, FOV_acquired, R, PF_Factor, radialBinWidthMode);
[regionA, phaseEncodeTable] = getRegionA( ...
    Matrix_Size_Acquired, FOV_acquired, pA, PF_Factor, R, radialBinWidthMode);

% Toggle this to "min" if you want shells quantized using the smaller
% phase/slice k-space pixel size.
disp("radialBinWidthMode = " + radialBinWidthMode)

radialBinGrid = zeros(Matrix_Size_Acquired(2), Matrix_Size_Acquired(3));
radialBinGrid(phaseEncodeTable.LinearIndex) = phaseEncodeTable.RadialBin;
retainedMask = false(Matrix_Size_Acquired(2), Matrix_Size_Acquired(3));
retainedMask(phaseEncodeTable.LinearIndex) = phaseEncodeTable.IsKeptAfterAcceleration;

figure('Name', 'TWIST Ordering Diagnostics', 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(regionA');
axis image;
set(gca, 'YDir', 'normal');
title('Region A');
xlabel('Phase Encode');
ylabel('Slice Encode');

nexttile;
imagesc(radialBinGrid');
axis image;
set(gca, 'YDir', 'normal');
title(sprintf('Radial Bins (%s dk)', radialBinWidthMode));
xlabel('Phase Encode');
ylabel('Slice Encode');
colorbar;

nexttile;
imagesc(retainedMask');
axis image;
set(gca, 'YDir', 'normal');
title('Predicted PF / GRAPPA Mask');
xlabel('Phase Encode');
ylabel('Slice Encode');

plotTWISTFrameOrdering( ...
    TWIST_sampling_order, Matrix_Size_Acquired, regionA, ...
    "TWIST B-Frame Ordering", 9);
