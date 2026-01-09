clc;
% Get current path
[pathstr,name,ext] = fileparts(mfilename('fullpath'));

% Add CIPG common package to path
disp('Adding mri breast virtual trial simulator package to MATLAB path...');
path(genpath([pathstr filesep() 'src']) ,path);
disp('Done')
