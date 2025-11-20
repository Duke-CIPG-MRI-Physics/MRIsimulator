function k_vec = computeKspaceAxis(FOV_mm, N)
% computeKspaceAxis  Compute 1D spatial frequency axis (cycles/mm)
%                    consistent with the DFT/FFT convention.
%
%   k_vec = computeKspaceAxis(FOV_mm, N)
%
%   Inputs:
%       FOV_mm : field of view in mm (scalar, > 0)
%       N      : number of samples along this dimension (integer, > 0)
%
%   Output:
%       k_vec  : 1 x N vector of spatial frequencies in cycles/mm,
%                ordered as fftshifted DFT frequencies, i.e.,
%                k = ((0:N-1) - floor(N/2)) / FOV_mm
%
%   Notes:
%       - This definition is consistent with using fftshift/ifftshift
%         around fftn/ifftn for image reconstruction.
%       - For even N, this is equivalent to (-N/2:N/2-1)/FOV_mm.

    arguments
        FOV_mm (1,1) double {mustBePositive}
        N      (1,1) double {mustBeInteger, mustBePositive}
    end

    % DFT-consistent frequencies (cycles/mm)
    k_vec = ((0:N-1) - floor(N/2)) / FOV_mm;
end
