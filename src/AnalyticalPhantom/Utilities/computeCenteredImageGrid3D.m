function [x_vec, y_vec, z_vec, X, Y, Z] = computeCenteredImageGrid3D(FOV_mm, N)
% computeCenteredImageGrid3D
%   Build 3D spatial grid (WORLD coordinates, mm), centered.
%
%   [x_vec, y_vec, z_vec, X, Y, Z] = computeCenteredImageGrid3D(FOV_mm, N)
%
%   Inputs:
%       FOV_mm : scalar or [1x3] double for [FOVx FOVy FOVz]
%       N      : [1x3] integer for [Nx Ny Nz]
%
%   Outputs:
%       x_vec,y_vec,z_vec : centered axes in mm
%       X,Y,Z             : ndgrid(x_vec,y_vec,z_vec)
%
%   Pairs with:
%       img = fftshift(ifftn(ifftshift(K)));

    arguments
        FOV_mm (1,:) double {mustBePositive}
        N      (1,3) double {mustBeInteger, mustBePositive}
    end

    % Allow scalar FOV→ replicate to 3 components
    if isscalar(FOV_mm)
        FOV_mm = repmat(FOV_mm, 1, 3);
    elseif numel(FOV_mm) ~= 3
        error('FOV_mm must be scalar or 1x3 [FOVx FOVy FOVz].');
    end

    Nx = N(1); Ny = N(2); Nz = N(3);

    x_vec = computeCenteredImageGrid1D(FOV_mm(1), Nx);
    y_vec = computeCenteredImageGrid1D(FOV_mm(2), Ny);
    z_vec = computeCenteredImageGrid1D(FOV_mm(3), Nz);

    [X, Y, Z] = ndgrid(x_vec, y_vec, z_vec);
end
