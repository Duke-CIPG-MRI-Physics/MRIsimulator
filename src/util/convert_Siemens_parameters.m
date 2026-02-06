function [FOV_acquired,kmatrix_size_complete,kmatrix_size_acquired,voxel_size_mm,...
    nyquist_resolution_mm,IMmatrix_crop_size] = convert_Siemens_parameters(Scan_Parameters)
%This function takes as input a Dictionary of Siemens scan parameters
%to output what the true dimensions and FOV of k-space are.


FOV_read_mm = Scan_Parameters("FOV read (mm)");
FOV_phase_pct = Scan_Parameters("FOV phase (%)");
oversampling_phase_pct = Scan_Parameters("Phase Oversampling (%)");
oversampling_slice_pct = Scan_Parameters("Slice Oversampling (%)");
slices_per_slab = Scan_Parameters("Slices per Slab");
slice_thickness_mm = Scan_Parameters("Slice Thickness (mm)");
base_resolution = Scan_Parameters("Base Resolution");
phase_resolution_pct = Scan_Parameters("Phase Resolution (%)");
slice_resolution_pct = Scan_Parameters("Slice Resolution (%)");



%% Derive oversampled FOV (before cropping)
phase_FOV_pct = FOV_phase_pct/100;
phase_oversampling = (1+(oversampling_phase_pct/100));
slice_oversampling = (1+(oversampling_slice_pct/100));

%this is the FOV truly captured
FOV_acquired = [FOV_read_mm, ...                                 % freq dir
    FOV_read_mm*phase_oversampling*phase_FOV_pct, ...            % phase dir
    slices_per_slab*slice_oversampling*slice_thickness_mm];      % slice dir

%% Derive oversampled matrix (before cropping, ignoring PI, PF, TWIST, etc)
% *NOTE* this also includes zero-padding from percent resolution,
 % this is the matrix we will ifft
kmatrix_size_complete = [base_resolution, ...              % freq dir
    base_resolution*phase_oversampling, ...          % phase dir
    slices_per_slab*slice_oversampling];             % slice dir


%% Check if matrix is integer
%  - sometimes oversampling is off from decimal rounding
if(~all(kmatrix_size_complete == ceil(kmatrix_size_complete)))
    phase_oversampling = ceil(base_resolution*phase_oversampling)/base_resolution;
    slice_oversampling = ceil(slices_per_slab*slice_oversampling)/slices_per_slab;

    oversampling_phase_pct = 100*(phase_oversampling-1);
    oversampling_slice_pct = 100*(slice_oversampling-1);

    kmatrix_size_complete = [base_resolution, ...     % freq dir
        base_resolution*phase_oversampling, ... % phase dir
        slices_per_slab*slice_oversampling];    % slice dir

end

%% Derive acquired matrix (before cropping, ignoring PI, PF, TWIST, etc)
% Note this excludes zero padding from percent phase/slice resolution
    %this is the matrix we will acquire
phase_resolution = (phase_resolution_pct/100);
slice_resolution = (slice_resolution_pct/100);

kmatrix_size_acquired = [kmatrix_size_complete(1), ...      % freq dir
    kmatrix_size_complete(2)*phase_resolution, ...    % phase dir
    kmatrix_size_complete(3)*slice_resolution];       % slice dir

% Acquired matrix must also be integer - somtimes phase/slice resolution
% makes this slightly off due to decimal reporting in pdf values
if(~all(kmatrix_size_acquired == ceil(kmatrix_size_acquired)))
    phase_resolution = ceil(kmatrix_size_complete(2)*phase_resolution)/kmatrix_size_complete(2);
    slice_resolution = ceil(kmatrix_size_complete(3)*slice_resolution)/kmatrix_size_complete(3);

    phase_resolution_pct = 100*phase_resolution;
    slice_resolution_pct = 100*slice_resolution;

    kmatrix_size_acquired = [kmatrix_size_complete(1), ...      % freq dir
        kmatrix_size_complete(2)*phase_resolution, ...    % phase dir
        kmatrix_size_complete(3)*slice_resolution];       % slice dir
    
end


% Derived voxel sizes (ignoring percent slice/phase resolution)
voxel_size_mm = FOV_acquired./kmatrix_size_complete;

% Derived nominal voxel sizes (considering percent slice/phase resolution)
nyquist_resolution_mm = FOV_acquired./kmatrix_size_acquired;

%% Derive Final Image Crop Size (Post-Reconstruction)
crop_read = base_resolution;
crop_phase = base_resolution * phase_FOV_pct;
crop_slice = slices_per_slab;

IMmatrix_crop_size = [crop_read, crop_phase, crop_slice];

% Sanity Check: Ensure dimensions are integers
if any(mod(IMmatrix_crop_size, 1) ~= 0)
    IMmatrix_crop_size = round(IMmatrix_crop_size);
end


%% Display basic parameters
% disp( '***Basic protocol parameters: ***');
% disp(['   FOV oversampled   (f x ph x sl):' num2str(FOV_oversampled(1)) ' x ' num2str(FOV_oversampled(2)) ' x ' num2str(FOV_oversampled(3)) '(mm)']);
% disp(['   Matrix acq os     (f x ph x sl):' num2str(matrix_acq_os(1)) ' x ' num2str(matrix_acq_os(2)) ' x ' num2str(matrix_acq_os(3))]);
% disp(['   Voxel size        (f x ph x sl):' num2str(voxel_size_mm(1)) ' x ' num2str(voxel_size_mm(2)) ' x ' num2str(voxel_size_mm(3)) '(mm)']);
% disp(['   Nyquist resolution (f x ph x sl):' num2str(nyquist_resolution_mm(1)) ' x ' num2str(nyquist_resolution_mm(2)) ' x ' num2str(nyquist_resolution_mm(3)) '(mm)']);



end