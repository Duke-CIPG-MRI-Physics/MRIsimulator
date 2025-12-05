function x_vec = computeCenteredImageGrid1D(FOV_mm, N)
% computeCenteredImageGrid1D
%   Return a centered spatial axis (WORLD coordinates, mm).
%
%   x_vec = computeCenteredImageGrid1D(FOV_mm, N)
%
%   Inputs:
%       FOV_mm : scalar double, field-of-view in mm
%       N      : scalar integer, number of samples
%
%   Output:
%       x_vec  : 1 x N vector of spatial coordinates [mm],
%                centered such that:
%                   x = 0 occurs at index floor(N/2)+1
%
%   Convention:
%       dx    = FOV_mm / N
%       x_vec = ((0:N-1) - floor(N/2)) * dx
%
%   This axis pairs with image data reconstructed via:
%       img = fftshift(ifftn(ifftshift(K)));

    arguments
        FOV_mm (1,1) double {mustBePositive}
        N      (1,1) double {mustBeInteger, mustBePositive}
    end

    dx = FOV_mm / N;
    x_vec = ((0:N-1) - floor(N/2)) * dx;
end
