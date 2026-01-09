function[downsampled_image] = FourierDownsample(input_image, fov, output_resolution)
% FOURIER_DOWNSAMPLE Downsamples an N-dimensional volume in the Fourier domain.
%   DOWNSAMPLED_IMAGE = FOURIER_DOWNSAMPLE(VOLUME, FOV, PIXEL_SIZE, OUTPUT_RESOLUTION)
%   reduces the spatial frequency content of an N-dimensional input image
%   volume by downsampling in the Fourier domain. The function retains the
%   original matrix size of the input volume but limits nonzero frequencies to
%   those corresponding to the desired resolution (aka it zero pads the
%   undersampled k-space. Note that this function will not change or rescale
%   the DC values, and since the output matrix size is the same as the
%   input matrix size, the intensities (signal) should not change on a
%   pixel-by pixel basis.
%
%   Inputs:
%       VOLUME            - N-dimensional input image volume.
%       FOV               - 1xN vector specifying the Field of View (FOV) in
%                           each dimension, e.g., [FOV_x, FOV_y, ..., FOV_n].
%       OUTPUT_RESOLUTION - 1xN vector specifying the desired output resolution
%                           in each dimension, e.g., [res_x, res_y, ..., res_n].
%
%   Outputs:
%       DOWNSAMPLED_IMAGE - N-dimensional image volume with downsampled
%                           frequencies, retaining the same matrix size as the
%                           input volume.
%
%   Example:
%       downsampled_image = Fourier_downsample(volume, [256, 256, 150], ...
%                                              [1, 1, 1.5], [2, 2, 3]);
%
%   See also FFTN, FFTSHIFT, IFFTN, IFFTSHIFT

% Get the original dimensions of the volume
volume_size = size(input_image);
num_dims = length(volume_size); % Determine the number of dimensions

% Calculate the number of frequencies needed in each dimension based on desired resolution
required_freqs = round(fov ./ output_resolution);
%
% % Calculate the number of frequencies needed in each dimension based on desired resolution
% required_freqs = round(fov ./ output_resolution);

% Calculate the center and the start/end indices for each dimension
start_idx = zeros(1, num_dims);
end_idx = zeros(1, num_dims);
for dim = 1:num_dims
    center = ceil(volume_size(dim) / 2); % Find the center index in each dimension
    start_idx(dim) = center - floor(required_freqs(dim) / 2); % Starting index for desired frequencies
    end_idx(dim) = start_idx(dim) + required_freqs(dim) - 1; % Ending index for desired frequencies
end

% Generate a cell array of indices for slicing in each dimension
slice_indices = arrayfun(@(dim) start_idx(dim):end_idx(dim), 1:num_dims, 'UniformOutput', false);

% Convert input volume to k-space and shift to center the frequencies
k_space_shifted = fftshift(fftn(input_image));

% Place the required frequencies into the center of the new k-space array
downsampled_k_space = zeros(volume_size);
downsampled_k_space(slice_indices{:}) = k_space_shifted(slice_indices{:});
% Calculate the center and the start/end indices for each dimension
start_idx = zeros(1, num_dims);
end_idx = zeros(1, num_dims);
for dim = 1:num_dims
    center = ceil(volume_size(dim) / 2); % Find the center index in each dimension
    start_idx(dim) = max(1,center - floor(required_freqs(dim) / 2)); % Starting index for desired frequencies
    end_idx(dim) = start_idx(dim) + required_freqs(dim) - 1; % Ending index for desired frequencies

    if(start_idx(dim) < 1)
        test = 1;
    end
    if(end_idx(dim) > volume_size(dim))
        test =2;
    end
end

% Generate a cell array of indices for slicing in each dimension
slice_indices = arrayfun(@(dim) start_idx(dim):end_idx(dim), 1:num_dims, 'UniformOutput', false);

if all(start_idx > 0) && all(end_idx <= volume_size)
    downsampled_k_space(slice_indices{:}) = k_space_shifted(slice_indices{:});
else
    % Debug information
    if(all(start_idx > 0))
        disp('All indices > 0');
    end
    if(all(end_idx <= volume_size))
        disp('All indices <= matrix size');
    end
    if(any(start_idx == 0))
        warning('zero indices exist!');
    end
    if(any(start_idx < 0))
        warning('Negative indices exist!');
    end
    if(any(end_idx == volume_size))
        warning('Index at matrix size!')
    end
    if(any(end_idx > volume_size))
        warning('Index exceeds matrix size!')
    end
    error('Invalid slice indices, please check FOV and output resolution values.');
end

% Place the required frequencies into the center of the new k-space array
downsampled_k_space = zeros(volume_size);
downsampled_k_space(slice_indices{:}) = k_space_shifted(slice_indices{:});

% Inverse shift and transform back to the spatial domain
downsampled_image = ifftn(ifftshift(downsampled_k_space));
end