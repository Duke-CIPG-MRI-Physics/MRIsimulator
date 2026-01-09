% Demo Explanation:
% This script demonstrates the generation of TWIST masks for a 3D k-space acquisition.
% The A region is sampled densely (white), while the B region is sampled sparsely in different frames (colored).
% Each frame is visualized individually with a unique color assigned to its B region.
% Finally, a combined RGB mask showing all frames is displayed.

FOV_m = [0.36 0.36 0.2];  % Field of view in meters for x, y, and z directions
matrixSize = [444 128 128];  % Matrix size in frequency, phase, and slice directions
A = 0.15;  % A region fraction for TWIST
B = 1/10;   % B region fraction for TWIST
TR_sec = 0.008;  % Repetition time in seconds
R = [3 2];  % Acceleration factors in phase and slice encoding directions
partialFourier = [7/8 6/8];  % Partial Fourier factors in phase and slice directions

% Generate TWIST masks and calculate frame and reconstruction times
[twistMasks, Tframe_sec, Trecon_sec] = calculateTwistMasks(FOV_m, matrixSize, A, B, TR_sec, R, partialFourier);

% Generate distinguishable colors for B frames, avoiding white and black
n_colors = max(twistMasks(:)) - 1;
colorV = distinguishable_colors(n_colors, [1 1 1; 0 0 0]);

% Create RGB mask with white for A region and different colors for B frames
rgbMask = zeros([size(twistMasks) 3]);
rgbMask(repmat(twistMasks == 1, [1 1 3])) = 1;  % Set A region to white

% Assign colors to B frames and add them to the RGB mask
for iC = 1:n_colors
    frameMask = twistMasks == (iC + 1);  % Mask for current B frame
    rgbMask(:, :, 1) = rgbMask(:, :, 1) + colorV(iC, 1) * frameMask;
    rgbMask(:, :, 2) = rgbMask(:, :, 2) + colorV(iC, 2) * frameMask;
    rgbMask(:, :, 3) = rgbMask(:, :, 3) + colorV(iC, 3) * frameMask;
end

% Determine subplot layout for visualizing individual frames
nRows = floor(sqrt(n_colors + 1));
nCols = ceil((n_colors + 1) / nRows);

% Visualize individual frames and the full RGB mask
figure();
for iC = 1:n_colors
    frameMask = zeros([size(twistMasks) 3]);
    frameMask(:, :, 1) = colorV(iC, 1) * (twistMasks == (iC + 1)) + (twistMasks == 1);
    frameMask(:, :, 2) = colorV(iC, 2) * (twistMasks == (iC + 1)) + (twistMasks == 1);
    frameMask(:, :, 3) = colorV(iC, 3) * (twistMasks == (iC + 1)) + (twistMasks == 1);

    subplot(nRows, nCols, iC);
    imshow(frameMask);
    axis image off;
    title(['Frame ' num2str(iC)]);
end

% Display the combined RGB mask for all frames
subplot(nRows, nCols, n_colors + 1);
imshow(rgbMask);
axis image off;
title('All Frames (Combined RGB Mask)');

figure()
imshow(rgbMask);
axis image off;
% title('All Frames (Combined RGB Mask)');
