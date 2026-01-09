clear;clc;close all;
% %% load sensitivities
% load('sens.mat')
% 
% %% generate single slice test image
% x = phantom(64);
% 
% %% multiply through sensitivity maps
% x = x.*sens;

load("xcat_MRI.mat")
x = abs(xcat_MRI_real_downsampled);
% replicate 3D output across 8 channels
% x = phantom(64)
% x = repmat(x, [1, 1, 64, 4]);



%% Fourier transform input

%2D
%y = fftshift(fft(fftshift(fft(x,[],1),1),[],2),2);
%3D
y = zeros(size(x));

for ii = 1:size(x,4)
    y(:,:,:,ii) = fftshift(fftn(x(:,:,:,ii)));
end

%% generate noisy slice aliased input, and calibration data
n_calib_rows = 24;
n_calib_columns = 24;
n_calib_slices = 24;

center_rows = round((size(y,1)-n_calib_rows)/2:((size(y,1)-n_calib_rows)/2)+n_calib_rows-1);
center_columns = round((size(y,2)-n_calib_columns)/2:((size(y,2)-n_calib_columns)/2)+n_calib_columns-1);
center_slices = round((size(y,3)-n_calib_slices)/2:((size(y,3)-n_calib_slices)/2)+n_calib_slices-1);

calib = permute(y,[4,1,2,3]);

calib = calib(:,center_rows,center_columns,center_slices,:);


%need to permute because grappa expect coils as first dim

input = permute(y,[4,1,2,3]);

%optionally add noise
%input = input + (randn(size(input)) + 1j*randn(size(input)))/10;

%Removing GRAPPA Lines
R = [1,2,2];
input(:,:,:,R(2):R(2):end) = 0;
input(:,:,R(3):R(3):end,:) = 0;


%% perform grappa reconstruction
out   = grappa(input, calib, R , [4,4,4]);

%undoing permutations
out = ipermute(out,[4,1,2,3]);
input = ipermute(input,[4,1,2,3]);


%% inverse Fourier transform and sum-of-squares combine to show outputs

for ii = 1:size(out,4)
     img_GRAPPA(:,:,:,ii) = ifftn(ifftshift(out(:,:,:,ii)));
     img_aliased(:,:,:,ii) = ifftn(ifftshift(input(:,:,:,ii)));
end

img_GRAPPA = sqrt(sum(abs(img_GRAPPA).^2, 4));
img_aliased = sqrt(sum(abs(img_aliased).^2, 4));


figure
imshow3D(squeeze(log(abs(input(1,:,:,:)))))
title('input')

figure
imshow3D(squeeze(log(abs(out(1,:,:,:)))))
title('out')

figure
imshow3D(abs(img_GRAPPA(:,:,:,1)))

figure
imshow3D(abs(img_aliased(:,:,:,1)))

