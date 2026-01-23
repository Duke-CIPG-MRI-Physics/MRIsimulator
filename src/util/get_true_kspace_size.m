function [FOV_oversampled,matrix_acq_os,voxel_size_mm,nominal_resolution_mm] = get_true_kspace_size(Scan_Parameters)
%This function takes as input a Dictionary of Siemens scan parameters
%to output what the true dimensions and FOV of k-space are.


% arguments (Input)
%     inputArg1
%     inputArg2
% end

FOV_read_mm = Scan_Parameters("FOV read (mm)");
FOV_phase_pct = Scan_Parameters("FOV phase (%)");
oversampling_phase_pct = Scan_Parameters("Phase Oversampling (%)");
oversampling_slice_pct = Scan_Parameters("Slice Oversampling (%)");
slices_per_slab = Scan_Parameters("Slices per Slab");
slice_thickness_mm = Scan_Parameters("Slice Thickness (mm)");
base_resolution = Scan_Parameters("Base Resolution (mm)");
phase_resolution_pct = Scan_Parameters("Phase Resolution (%)");
slice_resolution_pct = Scan_Parameters("Slice Resolution (%)");



% Derived oversampled FOV (before cropping)
FOV_phase_dec = FOV_phase_pct/100;
os_pe_dec = (1+(oversampling_phase_pct/100));
os_sl_dec = (1+(oversampling_slice_pct/100));
FOV_oversampled = [FOV_read_mm, ...                     % freq dir
    FOV_read_mm*os_pe_dec*FOV_phase_dec, ...            % phase dir
    slices_per_slab*os_sl_dec*slice_thickness_mm];      % slice dir

% Derived oversampled matrix (before cropping, ignoring PI, PF, TWIST, etc)
% *NOTE* this also includes zero-padding from percent phase/slice
% resolution
matrix_full_os = [base_resolution, ...      % freq dir
    base_resolution*os_pe_dec, ...
    slices_per_slab*os_sl_dec];

% Check if matrix is integer - sometimes oversampling is off from decimal
% rounding
if(~all(matrix_full_os == ceil(matrix_full_os)))
    os_pe_dec = ceil(base_resolution*os_pe_dec)/base_resolution;
    os_sl_dec = ceil(slices_per_slab*os_sl_dec)/slices_per_slab;

    oversampling_phase_pct = 100*(os_pe_dec-1);
    oversampling_slice_pct = 100*(os_sl_dec-1);

    matrix_full_os = [base_resolution, ...      % freq dir
        base_resolution*os_pe_dec, ...
        slices_per_slab*os_sl_dec];

  
end

% Derive acquired matrix (before cropping, ignoring PI, PF, TWIST, etc)
% Note this excludes zero padding from percent phase/slice resolution
pe_res_dec = (phase_resolution_pct/100);
sl_res_dec = (slice_resolution_pct/100);
matrix_acq_os = [matrix_full_os(1), ...
    matrix_full_os(2)*pe_res_dec, ...
    matrix_full_os(3)*sl_res_dec];

% Acquired matrix must also be integer - somtimes phase/slice resolution
% makes this slightly off due to decimal reporting in pdf values
if(~all(matrix_acq_os == ceil(matrix_acq_os)))
    pe_res_dec = ceil(matrix_full_os(2)*pe_res_dec)/matrix_full_os(2);
    sl_res_dec = ceil(matrix_full_os(3)*sl_res_dec)/matrix_full_os(3);

    phase_resolution_pct = 100*pe_res_dec;
    slice_resolution_pct = 100*sl_res_dec;

    matrix_acq_os = [matrix_full_os(1), ...
        matrix_full_os(2)*pe_res_dec, ...
        matrix_full_os(3)*sl_res_dec];


end


% Derived voxel sizes (ignoring percent slice/phase resolution)
voxel_size_mm = FOV_oversampled./matrix_full_os;

% Derived nominal voxel sizes (considering percent slice/phase resolution)
nominal_resolution_mm = FOV_oversampled./matrix_acq_os;

%% Display basic parameters
% disp( '***Basic protocol parameters: ***');
% disp(['   FOV oversampled   (f x ph x sl):' num2str(FOV_oversampled(1)) ' x ' num2str(FOV_oversampled(2)) ' x ' num2str(FOV_oversampled(3)) '(mm)']);
% disp(['   Matrix acq os     (f x ph x sl):' num2str(matrix_acq_os(1)) ' x ' num2str(matrix_acq_os(2)) ' x ' num2str(matrix_acq_os(3))]);
% disp(['   Voxel size        (f x ph x sl):' num2str(voxel_size_mm(1)) ' x ' num2str(voxel_size_mm(2)) ' x ' num2str(voxel_size_mm(3)) '(mm)']);
% disp(['   Nominal resolution (f x ph x sl):' num2str(nominal_resolution_mm(1)) ' x ' num2str(nominal_resolution_mm(2)) ' x ' num2str(nominal_resolution_mm(3)) '(mm)']);



end