clc;clear;close all

%% Input Data
%theFileName = 'C:\Users\Roberto\My Drive\Grad School\Research\RobbieTWIST\No Partial Fourier, Minimized GRAPPA\meas_MID00988_FID84788_t1_twist_R3_R2_NOPF_A4_B10.dat';

FileName = 'meas_MID00033_FID24847__UF_PF_A3_a30_b0.dat';

% Read the TWIX data
twix_objs = mapVBVD(FileName);
if(~iscell(twix_objs))
    twix_objs = {twix_objs};
end

% Get the k-space data
% 
% A_percentage = 

% {'Col'}  {'Cha'}    {'Lin'}    {'Par'}    {'Sli'}    {'Ave'}    {'Phs'}    {'Eco'}    {'Rep'}
%     {'Set'}    {'Seg'}    {'Ida'}    {'Idb'}    {'Idc'}    {'Idd'}    {'Ide'}


kspace_all = squeeze(twix_objs{2}.image(1,1,:,:,1,1,1,1,:,1,1,1,1,1,1,1));
kspace_all_binary = (kspace_all~=0);



imslice(kspace_all_binary);
%Code below is parameter specific


% kspace_filled = kspace_all_binary(:,:,1);
% 
% %If GRAPPA was used can can use the next line to remove blank rows
% kspace_grappa_binary = kspace_all_binary(3:3:end,2:2:end,:);
% 
% 
% %This is where we set the number of measurements based on B_percentage
% %Change the variable between grappa and all depending on visualization
% %needs
% kspace_full_sequence = kspace_all_binary(:,:,1:num_meas+1);
% 
% 
% figure
% imshow3D(kspace_full_sequence)
% title('kspace, full sequence, boolean')
% 
% A_region = uint8(sum(kspace_full_sequence,3)==num_meas);
% 
% kspace_full_sequence_double = uint8(kspace_full_sequence);
% B_region = uint8(zeros(size(kspace_full_sequence,1),size(kspace_full_sequence,2)));
% 
% for ii = 1:size(kspace_full_sequence,3)
% B_region = B_region + kspace_full_sequence_double(:,:,ii)*ii;
% end
% 
% B_region = B_region - (A_region*255);
% 
% 
% 
% %% 
% 
% figure
% tiledlayout('Horizontal')
% 
% nexttile
% imshow(A_region,[])
% 
% nexttile
% imagesc(B_region)
% colormap('copper')
% 
% true_A_percentage = sum(A_region(:))/numel(A_region)*100;
% 
% fprintf('True A%% = %f%%\n',true_A_percentage)
% 
% 
