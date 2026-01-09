%Mask Checking function, just opens stuff
clear;clc;close all;
load('Masks_V2.mat')

subplot(2,2,1)
imshow(fat_mask(:,:,100))
subplot(2,2,2)
imshow(fibroglandular_mask(:,:,100))
subplot(2,2,3)
imshow(lesion_mask(:,:,100))
subplot(2,2,4)
imshow(whole_breast_mask(:,:,100))

figure
sliceViewer(lesion_mask)
figure
sliceViewer(fibroglandular_mask)
figure
sliceViewer(boolean(fibroglandular_mask-lesion_mask))