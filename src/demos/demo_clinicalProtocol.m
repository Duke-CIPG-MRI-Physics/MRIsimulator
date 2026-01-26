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
base_resolution = 224; % 224 default
phase_resolution_pct = 100;
slice_resolution_pct = 80;
freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I (use 2 1 3 for R/L PE, S/I slice)

% FOV_read_mm = 350;
% FOV_phase_pct = 100;
% oversampling_phase_pct = 20; 
% oversampling_slice_pct = 33.3; 
% slices_per_slab = 240;
% slice_thickness_mm = 1; % Note this is the reconstructed slice thickness, not the nominal slice thickness 
% base_resolution = 224; % 224 default
% phase_resolution_pct = 100;
% slice_resolution_pct = 80;
% freq_phase_slice = [2 1 3]; % 1 = R/L, 2=A/P, 3 = S/I (use 2 1 3 for R/L PE, S/I slice)

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
FOV_oversampled_fps = [FOV_read_mm, ...                     % freq dir
    FOV_read_mm*os_pe_dec*FOV_phase_dec, ...            % phase dir
    slices_per_slab*os_sl_dec*slice_thickness_mm];      % slice dir

% Derived oversampled matrix (before cropping, ignoring PI, PF, TWIST, etc)
% *NOTE* this also includes zero-padding from percent phase/slice
% resolution
matrix_full_os_fps = [base_resolution, ...      % freq dir
    base_resolution*os_pe_dec, ...
    slices_per_slab*os_sl_dec];

% Check if matrix is integer - sometimes oversampling is off from decimal
% rounding
if(~all(matrix_full_os_fps == ceil(matrix_full_os_fps)))
    os_pe_dec = ceil(base_resolution*os_pe_dec)/base_resolution;
    os_sl_dec = ceil(slices_per_slab*os_sl_dec)/slices_per_slab;

    oversampling_phase_pct = 100*(os_pe_dec-1);
    oversampling_slice_pct = 100*(os_sl_dec-1);

    matrix_full_os_fps = [base_resolution, ...      % freq dir
        base_resolution*os_pe_dec, ...
        slices_per_slab*os_sl_dec];

    disp('Adjusting oversampling to give integer matrix...');
    disp(['   Adjusted oversampling_phase_pct=' num2str(oversampling_phase_pct) '%'])
    disp(['   Adjusted oversampling_slice_pct=' num2str(oversampling_slice_pct) '%'])
end

% Derive acquired matrix (before cropping, ignoring PI, PF, TWIST, etc)
% Note this excludes zero padding from percent phase/slice resolution
pe_res_dec = (phase_resolution_pct/100);
sl_res_dec = (slice_resolution_pct/100);
matrix_acq_os_fps = [matrix_full_os_fps(1), ...
    matrix_full_os_fps(2)*pe_res_dec, ...
    matrix_full_os_fps(3)*sl_res_dec];

% Acquired matrix must also be integer - somtimes phase/slice resolution
% makes this slightly off due to decimal reporting in pdf values
if(~all(matrix_acq_os_fps == ceil(matrix_acq_os_fps)))
    pe_res_dec = ceil(matrix_full_os_fps(2)*pe_res_dec)/matrix_full_os_fps(2);
    sl_res_dec = ceil(matrix_full_os_fps(3)*sl_res_dec)/matrix_full_os_fps(3);

    phase_resolution_pct = 100*pe_res_dec;
    slice_resolution_pct = 100*sl_res_dec;

    matrix_acq_os_fps = [matrix_full_os_fps(1), ...
        matrix_full_os_fps(2)*pe_res_dec, ...
        matrix_full_os_fps(3)*sl_res_dec];

    disp('Adjusting phase/slice resolution to give integer matrix...');
    disp(['   Adjusted phase_resolution_pct=' num2str(phase_resolution_pct) '%'])
    disp(['   Adjusted slice_resolution_pct=' num2str(slice_resolution_pct) '%'])
end


% Derived voxel sizes (ignoring percent slice/phase resolution)
voxel_size_mm = FOV_oversampled_fps./matrix_full_os_fps;

% Derived nominal voxel sizes (considering percent slice/phase resolution)
nominal_resolution_mm = FOV_oversampled_fps./matrix_acq_os_fps;

%% Display basic parameters
disp( '***Basic protocol parameters: ***');
disp(['   FOV oversampled   (f x ph x sl):' num2str(FOV_oversampled_fps(1)) ' x ' num2str(FOV_oversampled_fps(2)) ' x ' num2str(FOV_oversampled_fps(3)) '(mm)']);
disp(['   Matrix acq os     (f x ph x sl):' num2str(matrix_acq_os_fps(1)) ' x ' num2str(matrix_acq_os_fps(2)) ' x ' num2str(matrix_acq_os_fps(3))]);
disp(['   Voxel size        (f x ph x sl):' num2str(voxel_size_mm(1)) ' x ' num2str(voxel_size_mm(2)) ' x ' num2str(voxel_size_mm(3)) '(mm)']);
disp(['   Nominal resoluton (f x ph x sl):' num2str(nominal_resolution_mm(1)) ' x ' num2str(nominal_resolution_mm(2)) ' x ' num2str(nominal_resolution_mm(3)) '(mm)']);

% Robbie - go back to TWIST here

%% Order k-space 
% * NOTE * this does not use TWIST, PI, PF yet so its much slower!
[t_s] = orderRectilinearKspace(matrix_acq_os_fps, dt_s, TR_s); % Scott - rename to calculateKspaceSampleTimes_rectilinear

encodingFullStr = formatEncodingString(freq_phase_slice);
disp(['   Encoding          (f x ph x sl):' encodingFullStr])

%% 5) Construct the breast phantom with the embedded enhancing vessel
disp('Constructing phantom');
breastPhantomParams = createBreastPhantomParams();
phantom = BreastPhantom(breastPhantomParams);

%% 6) Visualize contrast plug-flow timing relative to acquisition
t_sorted_s = sort(t_s(:));
contrastLength_mm = calculatePlugFlowInVessels(t_sorted_s, breastPhantomParams);
acqStart_s = min(t_sorted_s);
acqEnd_s = max(t_sorted_s);
rectTime_s = [acqStart_s acqStart_s acqEnd_s acqEnd_s];
rectAmp = [0 1 1 0];

figure('Name', 'Contrast Plug Flow and Acquisition Window');
yyaxis left
plot(t_sorted_s, contrastLength_mm, 'LineWidth', 1.5);
ylabel('Contrast length [mm]')
xlabel('Time [s]')

yyaxis right
plot(rectTime_s, rectAmp, 'LineWidth', 1.5);
ylabel('K-space acquisition window')
ylim([-0.1 1.1])
grid on


%% 3) Build WORLD k-space grid and map to the rectilinear ordering
disp('Constructing k-space grid');
[~, ~, ~, kfreq, kPhase, kSlice] = computeKspaceGrid3D(FOV_oversampled_fps, matrix_acq_os_fps);
k_fps = [kfreq(:) kPhase(:) kSlice(:)]';

%% Permute dimmensions to convert FPS to XYZ
disp('Permuting');
kspaceSize = size(kfreq);
[k_xyz, fps_to_xyz] = mapKspaceFpsToXyz(k_fps, freq_phase_slice);
clear k_fps kfreq kPhase kSlice;

%% 6) Compute analytic k-space for the phantom in ordered acquisition space
fprintf('Evaluating analytic k-space...\n');
K = phantom.kspaceAtTime(k_xyz(1,:), k_xyz(2,:), k_xyz(3,:), t_s);
clear t_s;
K = reshape(K,matrix_acq_os_fps);

%% Permute so that image is always oriented with dims:
% 1 = R/L (x)
% 2 = A/P (y)
% 3 = S/II (z)
K = permute(K, fps_to_xyz); 

%% 6) Reconstruct 3D image via inverse FFT
fprintf('Performing 3D inverse FFT...\n');
clear k_xyz phantom;
img_viaKspace = fftshift(ifftn(ifftshift(K)));

%% display
imslice(abs(img_viaKspace))
