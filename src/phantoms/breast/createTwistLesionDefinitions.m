function lesions = createTwistLesionDefinitions(intensityFunction)
% createTwistLesionDefinitions  Fixed TWIST lesion layout in right-breast mm coordinates.
%   lesions = createTwistLesionDefinitions(intensityFunction) returns the
%   four lesion definitions used by the analytical TWIST demo and
%   simulator. Lesion centers are stored as explicit millimeter shifts
%   relative to the center of the right breast.
%
%   The fixed millimeter shifts were derived once from the historical
%   voxel-scaled lesion centers using the clinical ultrafast reference
%   protocol in Breast_Ultrafast_scan_parameters.mat:
%       - reference voxel size [freq, phase, slice] = [2.1875, 2.1875, 1] mm
%       - the original demo applied those voxel counts directly to
%         center_mm [x y z], so the preserved fixed shifts below follow
%         that same world-coordinate convention
%
%   Inputs
%   ------
%   intensityFunction : function_handle, optional
%       Lesion intensity function handle evaluated by BreastPhantom.
%
%   Output
%   ------
%   lesions : [4x1] struct
%       Fields:
%           label              - lesion size label: L, M, S, XS
%           mm_shift           - [1x3] fixed millimeter shift [x y z]
%           center_mm          - alias of mm_shift for BreastPhantom
%           radius_mm          - spherical lesion radius [mm]
%           intensityFunction  - intensity function handle

if nargin < 1 || isempty(intensityFunction)
    params = createBreastPhantomParams();
    intensityFunction = params.lesionIntensityFunction;
elseif ~isa(intensityFunction, 'function_handle')
    error('createTwistLesionDefinitions:InvalidIntensityFunction', ...
        'intensityFunction must be a function handle.');
end

lesionLabels = {'L'; 'M'; 'S'; 'XS'};
lesionShifts_mm = [ ...
    -40, 0, 0; ...
    39.375, 15.625, 16; ...
    -21.875, -21.875, -14; ...
    13.125, -43.75, 22];
lesionRadii_mm = [10; 5; 2.5; 1.25];

nLesions = numel(lesionRadii_mm);
lesions = repmat(struct( ...
    'label', '', ...
    'mm_shift', zeros(1, 3), ...
    'center_mm', zeros(1, 3), ...
    'radius_mm', 0, ...
    'intensityFunction', intensityFunction), nLesions, 1);

for lesionIdx = 1:nLesions
    lesions(lesionIdx).label = lesionLabels{lesionIdx};
    lesions(lesionIdx).mm_shift = lesionShifts_mm(lesionIdx, :);
    lesions(lesionIdx).center_mm = lesionShifts_mm(lesionIdx, :);
    lesions(lesionIdx).radius_mm = lesionRadii_mm(lesionIdx);
    lesions(lesionIdx).intensityFunction = intensityFunction;
end
end
