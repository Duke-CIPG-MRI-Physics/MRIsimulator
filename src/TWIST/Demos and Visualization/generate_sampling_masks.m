function [kspace_mask, calibration_mask] = generate_sampling_masks(matrix_size, R, partial_fourier_factors, num_calibration_lines,showFigure)
    % Generate masks for 3D k-space sampling with parallel imaging and partial Fourier
    %
    % Inputs:
    % - matrix_size: [Nx, Ny, Nz] (frequency, phase, slice)
    % - R: [R_phase, R_slice] acceleration factors for phase and slice encodes
    % - partial_fourier_factors: [PF_phase, PF_slice] partial Fourier factors for phase and slice
    % - num_calibration_lines: [Cal_phase, Cal_slice] number of calibration lines in phase and slice directions
    %
    % Outputs:
    % - kspace_mask: 3D k-space sampling mask
    % - calibration_mask: 3D mask for the calibration region only
    
    % Unpack matrix size and parameters
    Nx = matrix_size(1);
    Ny = matrix_size(2); % Phase encodes
    Nz = matrix_size(3); % Slice encodes
    R_phase = R(1);
    R_slice = R(2);
    PF_phase = partial_fourier_factors(1);
    PF_slice = partial_fourier_factors(2);
    Cal_phase = num_calibration_lines(1);
    Cal_slice = num_calibration_lines(2);
    
    % Initialize masks
    kspace_mask = zeros(Ny, Nz);
    calibration_mask = zeros(Ny, Nz);
    
    % Compute number of partial Fourier lines for phase and slice
    num_partial_phase = round(Ny * PF_phase);
    num_partial_slice = round(Nz * PF_slice);
    
    % Compute center regions for calibration lines
    center_phase_start = floor(Ny / 2) - floor(Cal_phase / 2) + 1;
    center_phase_end = center_phase_start + Cal_phase - 1;
    center_slice_start = floor(Nz / 2) - floor(Cal_slice / 2) + 1;
    center_slice_end = center_slice_start + Cal_slice - 1;
    
    % Apply R acceleration with partial Fourier sampling
    kspace_mask(1:R_phase:num_partial_phase, 1:R_slice:num_partial_slice) = 1;
    
    % Add calibration lines in the center
    % kspace_mask(center_phase_start:center_phase_end, center_slice_start:center_slice_end) = 1;
    calibration_mask(center_phase_start:center_phase_end, center_slice_start:center_slice_end) = 1;
    
    if(showFigure)
        % Display the masks
        figure;
        subplot(1, 2, 1);
        imagesc(kspace_mask);
        colormap(gray);
        axis equal tight;
        title('K-space Sampling Mask (Phase vs Slice)');

        subplot(1, 2, 2);
        imagesc(calibration_mask);
        colormap(gray);
        axis equal tight;
        title('Calibration Phase vs Slice Mask');
    end
end
