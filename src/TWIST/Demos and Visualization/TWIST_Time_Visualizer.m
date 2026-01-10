clear;clc;close all;

pA = .01:.01:.5;
N = 2:1:100;

PE_plane_size = [224,224];
TR = 5.88e-3;

% Create grids for all combinations
[N_grid, pA_grid] = meshgrid(N, pA);

% Vectorized calculation (note the '.' before operators for element-wise ops)
frame_time = pA_grid.^2 + ((1 - pA_grid.^2) ./ N_grid);

% Plot
surf(N_grid, pA_grid, frame_time);
xlabel('N')
ylabel('pA')
zlabel('temporal resolution')

% 1. Calculate the Gradient
% gradient(F, hx, hy) takes the spacing of your x (N) and y (pA) axes
% hx = spacing between N values (1)
% hy = spacing between pA values (0.01)
spacing_N = N(2) - N(1);   
spacing_pA = pA(2) - pA(1); 

[dTime_dN, dTime_dpA] = gradient(frame_time, spacing_N, spacing_pA);

% 2. Calculate Magnitude (Steepness)
% This combines the change from both N and pA into one "slope" value
slope_magnitude = sqrt(dTime_dN.^2 + dTime_dpA.^2);

% 3. Plot
figure('Name', 'Slope Analysis');
surf(N, pA, slope_magnitude);
xlabel('N (Undersampling)');
ylabel('pA (Central Region)');
zlabel('Gradient Magnitude (s)');
title('Steepness of Time Change (Gradient Magnitude)');
colorbar; % Adds a color scale for easier reading


