function k_vec = computeKspaceGrid1D(FOV_mm, N)
% computeKspaceGrid1D  1D k-space axis (WORLD coords, cycles/mm).
%
%   k_vec = computeKspaceGrid1D(FOV_mm, N)
%
%   Inputs:
%       FOV_mm : scalar double, field-of-view in mm
%       N      : scalar integer, number of samples
%
%   Output:
%       k_vec  : 1 x N vector of spatial frequencies [cycles/mm],
%                centered such that k=0 is in the middle index after
%                fftshift/ifftshift. This pairs with:
%                   img = fftshift(ifftn(ifftshift(K))).
%
%   Convention:
%       k = ((0:N-1) - floor(N/2)) / FOV_mm
%
%   So Δk = 1/FOV_mm and the total k-range is ~[-N/(2*FOV), +N/(2*FOV)).

    arguments
        FOV_mm (1,1) double {mustBePositive}
        N      (1,1) double {mustBeInteger, mustBePositive}
    end

    idx   = (0:N-1) - floor(N/2);  % centered integer indices
    k_vec = idx / FOV_mm;          % cycles/mm
end
