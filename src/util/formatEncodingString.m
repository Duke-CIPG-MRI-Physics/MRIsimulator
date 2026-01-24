function encodingStr = formatEncodingString(Freq_Phase_Slice_dirs)
% formatEncodingString Format encoding direction labels for display.
%
%   encodingStr = formatEncodingString(Freq_Phase_Slice_dirs)
%
%   Inputs:
%       Freq_Phase_Slice_dirs : [1x3] integer, encoding directions for
%                        [frequency, phase, slice].
%                        1 = R/L, 2 = A/P, 3 = S/I.
%
%   Output:
%       encodingStr : char row vector of formatted encoding directions.
%
%   Example:
%       encodingStr = formatEncodingString([2 1 3]);

    arguments
        Freq_Phase_Slice_dirs (1,3) double {mustBeInteger, mustBePositive}
    end

    validateattributes(Freq_Phase_Slice_dirs, {'double'}, ...
        {'integer', '>=', 1, '<=', 3, 'numel', 3}, mfilename, 'Freq_Phase_Slice_dirs');

    labels = {'R/L', 'A/P', 'S/I'};

    label_values = labels(Freq_Phase_Slice_dirs);
    encodingStr = sprintf('freq: %s, phase: %s, slice: %s', ...
        label_values{:});
end
