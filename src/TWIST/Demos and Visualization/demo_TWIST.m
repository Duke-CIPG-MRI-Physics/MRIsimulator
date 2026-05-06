clear; clc; close all;
% Dedicated TWIST ordering diagnostic demo.

% --- Inputs ---
pA = 0.04; 
pB = 0.1; 
Matrix_Size_Acquired = [1,96,96];
FOV_acquired = [1,350,96];
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


%%

% --- Correct Visualization Setup ---
nRows = Matrix_Size_Acquired(2);
nCols = Matrix_Size_Acquired(3);
nFrames = max(TWIST_sampling_order.Frame) + 1; % +1 because frames start at 0

% Initialize 3D mask: [Phase, Slice, Time]
sampled_mask = false(nRows, nCols, nFrames);

% Loop through each row of the table to place points in the correct "Time" slice
for i = 1:height(TWIST_sampling_order)
    r = TWIST_sampling_order.("Row (phase)")(i);
    c = TWIST_sampling_order.("Column (slice)")(i);
    f = TWIST_sampling_order.Frame(i) + 1; % Shift 0-index to 1-index for MATLAB

    sampled_mask(r, c, f) = true;
end

sampled_mask_swirl = zeros(size(sampled_mask(:,:,1)));

for i =2:size(sampled_mask,3)
    sampled_mask_swirl = sampled_mask_swirl + (double(sampled_mask(:,:,i)) * i);
end
sampled_mask_swirl(sampled_mask_swirl == max(sampled_mask_swirl)) = 0;
figure
imshow(sampled_mask_swirl,[])

sliceViewer(double(sampled_mask)); 
title('TWIST Sampling Mask (Scroll through Frames)');