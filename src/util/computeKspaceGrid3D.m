function [kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N)
% computeKspaceGrid3D  3D k-space grid (WORLD coords, cycles/mm).
%
%   [kx_vec, ky_vec, kz_vec, kx, ky, kz] = computeKspaceGrid3D(FOV_mm, N)
%
%   Inputs:
%       FOV_mm : scalar or [1x3] double, FOV in mm for [x y z].
%                 - If scalar, the same FOV is used for all three axes.
%       N      : [1x3] integer, [Nx Ny Nz] matrix size.
%
%   Outputs:
%       kx_vec, ky_vec, kz_vec : 1D frequency axes [cycles/mm]
%       kx, ky, kz             : 3D ndgrid arrays of k-space coordinates.
%
%   Designed to pair with:
%       img = fftshift(ifftn(ifftshift(K)));

    arguments
        FOV_mm (1,:) double {mustBePositive}
        N      (1,3) double  {mustBeInteger, mustBePositive}
    end

    % Allow scalar FOV → same in all directions
    if isscalar(FOV_mm)
        FOV_mm = repmat(FOV_mm, 1, 3);
    elseif numel(FOV_mm) ~= 3
        error('FOV_mm must be scalar or 1x3 [FOVx FOVy FOVz] in mm.');
    end

    kx_vec = computeKspaceGrid1D(FOV_mm(1), N(1));
    ky_vec = computeKspaceGrid1D(FOV_mm(2), N(2));
    kz_vec = computeKspaceGrid1D(FOV_mm(3), N(3));

    % Note that matlab uses (x,y) notation for meshgrid/ndgrid, but indexes 
    % according to (row,col), so dims 1 and two have to be swapped to 
    % ensure sampling is in the correct order
    [kx, ky, kz] = ndgrid(kx_vec, ky_vec, kz_vec);
end
