function [kx, ky, kz] = calculate_spatial_frequencies(FOV_m, matrixSize)
    % Calculate spatial frequencies for an fftshifted fftn of an image
    %
    % Inputs:
    % - FOV: [FOVx, FOVy, FOVz] (2D or 3D vector specifying field of view in each direction)
    % - matrixSize: [Nx, Ny, Nz] (2D or 3D vector specifying matrix size in each direction)
    %
    % Outputs:
    % - kx, ky, kz: 1D arrays representing spatial frequencies (in units of inverse FOV)
    %               kz is empty for 2D inputs

    % Determine dimensionality (2D or 3D)
    dim = length(FOV_m);

    % Calculate spatial frequency axes
    kx = generate_frequency_axis(FOV_m(1), matrixSize(1));
    ky = generate_frequency_axis(FOV_m(2), matrixSize(2));
    
    if dim == 3
        kz = generate_frequency_axis(FOV_m(3), matrixSize(3));
    else
        kz = [];
    end
end

function k = generate_frequency_axis(FOV, N)
    % Helper function to generate a frequency axis for a given FOV and matrix size
    %
    % Inputs:
    % - FOV: Field of view in a specific dimension
    % - N: Matrix size in that dimension
    %
    % Output:
    % - k: 1D array of spatial frequencies

    dk = 1 / FOV;  % Frequency step size
    if mod(N, 2) == 0
        % Even-sized matrix: DC at center
        k = (-N/2:N/2-1) * dk;
    else
        % Odd-sized matrix: DC at center
        k = (-(N-1)/2:(N-1)/2) * dk;
    end
end
