clear all; 
close all; 
clc; 
% restoredefaultpath; 
% run('C:\code\mriSimulator\cipg_setup.m');
% savepath

% parameters from clinical ultrafast protocol
FOV_read_mm = 350;
FOV_phase_pct = 100;
oversampling_phase_pct = 20; 
oversampling_slice_pct = 33.3; 
slices_per_slab = 240;
slice_thickness_mm = 1; % Note this is the reconstructed slice thickness, not the nominal slice thickness 
base_resolution = 224;
phase_resolution_pct = 100;
slice_resolution_pct = 80;
freq_phase_slice = [3 2 1]; % 1 = R/L, 2=A/P, 3 = S/I

% TODO - interpolation is off in current protocol, but we could consider it as an option in the future...

% Contrast parameters
rBW_HzPerPix = 570;
TR_s = 5.88E-3;   
TE_s = 2.63E-3;

% Derived contrast paramters
rBW_Hz = rBW_HzPerPix*base_resolution;
dt_s = 1/rBW_Hz;   % dwell time between frequency-encode samples [s]

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

    disp('Adjusting oversampling to give integer matrix...');
    disp(['   Adjusted oversampling_phase_pct=' num2str(oversampling_phase_pct) '%'])
    disp(['   Adjusted oversampling_slice_pct=' num2str(oversampling_slice_pct) '%'])
end

% Derived acquired matrix (before cropping, ignoring PI, PF, TWIST, etc)
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

    disp('Adjusting phase/slice resolution to give integer matrix...');
    disp(['   Adjusted phase_resolution_pct=' num2str(phase_resolution_pct) '%'])
    disp(['   Adjusted slice_resolution_pct=' num2str(slice_resolution_pct) '%'])
end


% Derived voxel sizes (ignoring percent slice/phase resolution)
voxel_size_mm = FOV_oversampled./matrix_full_os;

% Derived nominal voxel sizes (considering percent slice/phase resolution)
nominal_resolution_mm = FOV_oversampled./matrix_acq_os;

%% Display basic parameters
disp( '***Basic protocol parameters: ***');
disp(['   FOV oversampled   (f x ph x sl):' num2str(FOV_oversampled(1)) ' x ' num2str(FOV_oversampled(2)) ' x ' num2str(FOV_oversampled(3)) '(mm)']);
disp(['   Matrix acq os     (f x ph x sl):' num2str(matrix_acq_os(1)) ' x ' num2str(matrix_acq_os(2)) ' x ' num2str(matrix_acq_os(3))]);
disp(['   Voxel size        (f x ph x sl):' num2str(voxel_size_mm(1)) ' x ' num2str(voxel_size_mm(2)) ' x ' num2str(voxel_size_mm(3)) '(mm)']);
disp(['   Nominal resoluton (f x ph x sl):' num2str(nominal_resolution_mm(1)) ' x ' num2str(nominal_resolution_mm(2)) ' x ' num2str(nominal_resolution_mm(3)) '(mm)']);


%% Order k-space 
% * NOTE * this does not use TWIST, PI, PF yet so its much slower!
[kOrderedIdx, t_s] = orderRectilinearKspace(matrix_acq_os, freq_phase_slice, dt_s, TR_s);
t_s = t_s(:); % use sampling timestamps as the phantom time base

encodingFullStr = '';
for iDim = 1:3
    switch freq_phase_slice(iDim)
        case 1
            thisDimStr = 'R/L';
        case 2
            thisDimStr = 'A/P';
        case 3
            thisDimStr = 'S/I';
        otherwise
            error('freq_phase_slice can only be numbers 1-3')
    end
    if(iDim == 1)
        encodingFullStr = thisDimStr;
    else
        encodingFullStr = [encodingFullStr ' x ' thisDimStr];
    end
end

disp(['   Encoding          (f x ph x sl):' encodingFullStr])

%% 5) Construct the breast phantom with the embedded enhancing vessel
phantom = BreastPhantom(t_s);
clear t_s;

%% 3) Build WORLD k-space grid and map to the rectilinear ordering
[kx_vec, ky_vec, kz_vec, ~, ~, ~] = computeKspaceGrid3D(FOV_oversampled, matrix_acq_os);
kx_orderedIdx = kOrderedIdx(1, :);
ky_orderedIdx = kOrderedIdx(2, :);
kz_orderedIdx = kOrderedIdx(3, :);
clear kOrderedIdx;

ordKsx_kx = kx_vec(kx_orderedIdx);
ordKsx_ky = ky_vec(ky_orderedIdx); 
ordKsx_kz = kz_vec(kz_orderedIdx);
clear kx_vec ky_vec kz_vec;

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');
K_ordered = phantom.kspace(ordKsx_kx, ordKsx_ky, ordKsx_kz);
clear ordKsx_kx ordKsx_ky ordKsx_kz phantom;

% Reassemble onto the kx/ky/kz grid for reconstruction
K = zeros(matrix_acq_os);
linearIdx = sub2ind(matrix_acq_os, kx_orderedIdx, ky_orderedIdx, kz_orderedIdx);
clear kx_orderedIdx ky_orderedIdx kz_orderedIdx;
K(linearIdx) = K_ordered;
clear K_ordered;

%% 6) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
img_viaKspace = fftshift(ifftn(ifftshift(K)));



%% display
imslice(abs(img_viaKspace))