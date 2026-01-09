%GRE Simulator for use with XCat Breast Phantoms
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

fg_signal = MTrans;

%Adipose Tissue
T1_s = 423e-3;
T2star_s = 36.36e-3;
dfreq_Hz = 440;

[~,Mss_TE,~] = GRE_steadyState(flip_angle_radians,T1_s,T2star_s,TE_s,TR_s,dfreq_Hz);

[MTrans,~] = calculateMxyAndMz(Mss_TE);

adp_signal = MTrans;



%% --- Importing Xcat

load("XCAT_CTA0296.mat");


background_mask = Breast_Phantom == 0;
adp_mask = Breast_Phantom == 1;
fg25_mask = Breast_Phantom == 2;
fg50_mask = Breast_Phantom == 3;
fg75_mask = Breast_Phantom == 4;
fg_mask = Breast_Phantom == 5;
skin_mask = Breast_Phantom == 6;

%Calculating weighted signal values
Percent_fg = statistics{2,2:7};

fg25_signal = Percent_fg(2)*fg_signal + (1-Percent_fg(2))*adp_signal;
fg50_signal = Percent_fg(3)*fg_signal + (1-Percent_fg(3))*adp_signal;
fg75_signal = Percent_fg(4)*fg_signal + (1-Percent_fg(4))*adp_signal;
skin_signal = Percent_fg(5)*fg_signal;

%Multiplying masks by respective signal values
ADP_image = adp_mask.*adp_signal;
fg25_image = fg25_mask.*fg25_signal;
fg50_image = fg50_mask.*fg50_signal;
fg75_image = fg75_mask.*fg75_signal;
FG_image = fg_mask.*fg_signal;
skin_image = skin_mask.*skin_signal;


%% --- Combining scaled masks
Combined_output = ADP_image + fg25_image + fg50_image + fg75_image + FG_image + skin_image;
Combined_output = single(Combined_output);

figure
imshow3D(Combined_output);


%% --- Trying to make it stupid 4D

%Gonna be really stupid and just make the pure fg_mask slowly get brighter,
%updating only for the amount of measurements, might be a good way to make
%function based in the future

Combined_output_4D = Combined_output;

Num_Measurements = 6;

fg_signal = linspace(0,.3,Num_Measurements+1);

xcat_MRI_4D = zeros([size(Combined_output_4D),Num_Measurements+1]);

for ii = 1:Num_Measurements+1
    FG_image = fg_mask.*fg_signal(ii);
    Combined_output_4D(fg_mask) = FG_image(fg_mask);
    xcat_MRI_4D(:,:,:,ii) = Combined_output_4D; 
end

xcat_MRI_4D = single(xcat_MRI_4D);

figure
imshow3D(squeeze(xcat_MRI_4D(:,:,200,:)))