%Testing out applying GRE_steadyState.m to masks
clear;clc;close all

%Sequence parameters
flip_angle_degrees = 20;
flip_angle_radians = flip_angle_degrees*pi/180;
TE_s = 2.46e-3;
TR_s = 6e-3;

%Fibroglandular Tissue
T1_s = 1680e-3;
T2star_s = 30e-3;
dfreq_Hz = 0;

[~,Mss_TE,~] = GRE_steadyState(flip_angle_radians,T1_s,T2star_s,TE_s,TR_s,dfreq_Hz);

[MTrans,~] = calculateMxyAndMz(Mss_TE);

FG_signal = MTrans;

%Adipose Tissue
T1_s = 423e-3;
T2star_s = 36.36e-3;
dfreq_Hz = 440;

[~,Mss_TE,~] = GRE_steadyState(flip_angle_radians,T1_s,T2star_s,TE_s,TR_s,dfreq_Hz);

[MTrans,~] = calculateMxyAndMz(Mss_TE);

ADP_signal = MTrans;

%Lesion
T1_s = 1000e-3;
T2star_s = 30e-3;
dfreq_Hz = 0;

[~,Mss_TE,~] = GRE_steadyState(flip_angle_radians,T1_s,T2star_s,TE_s,TR_s,dfreq_Hz);

[MTrans,~] = calculateMxyAndMz(Mss_TE);

Lesion_signal = MTrans;

%Combining Values with masks
tic
load('Masks_V2.mat')
toc

%Converting boolean masks to singles
fibroglandular_mask = single(fibroglandular_mask);
fat_mask = single(fat_mask);
lesion_mask = single(lesion_mask);

%Multiplying masks by respective signal values
FG_image = fibroglandular_mask.*FG_signal;
ADP_image = fat_mask.*ADP_signal;
LE_image = lesion_mask.*Lesion_signal;

%Combining scaled masks
Combined_output = FG_image + ADP_image + LE_image;

sliceViewer(Combined_output)
