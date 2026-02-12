%This script is designed for performance testing on simulated
%Ultrafast breast MRI

clear;clc;close all;
load('Phantom_Output.mat')

%% Plotting 

%Ground truth
gt_timing = .01:.01:max(phantom_simulated.outputs.timing); %ground truth timing
gt_at_TWIST_pts = phantom_simulated.inputs.breastPhantomParams.lesionIntensityFunction(phantom_simulated.outputs.timing)...
    +phantom_simulated.inputs.breastPhantomParams.breastIntensity;


figure;
plot(gt_timing,phantom_simulated.inputs.breastPhantomParams.lesionIntensityFunction(gt_timing) ...
    +phantom_simulated.inputs.breastPhantomParams.breastIntensity);

%TODO: build function to output lesion ROI
x_range = 65:71;
y_range = 68:74;
z_range = 157:164;


roi_volume = phantom_simulated.outputs.phantom_magnitude(x_range,y_range,z_range,:);
roi_mean = mean(roi_volume, [1 2 3]);
contrast_values_measured = squeeze(roi_mean);

hold on
plot(phantom_simulated.outputs.timing,abs(contrast_values_measured),'.-','MarkerSize',15)
scatter(phantom_simulated.outputs.timing,gt_at_TWIST_pts,...
    15,'MarkerFaceColor',[1, 0.5, 0],'MarkerEdgeColor',[1, 0.5, 0])
legend("Ground Truth","TWIST Measured")
hold off

title("Contrast Wash-in")
xlabel("Time (s)")
ylabel("Pixel Value")

%% - Deviation

percent_deviation_from_GT = 100*abs((contrast_values_measured-gt_at_TWIST_pts)./gt_at_TWIST_pts);

max_percent_deviation_from_GT = max(percent_deviation_from_GT);
fprintf("Max deviation from ground-truth: %g%%\n",max_percent_deviation_from_GT)
average_percent_deviation_from_GT = mean(percent_deviation_from_GT);
fprintf("Avg. deviation from ground-truth: %g%%\n",average_percent_deviation_from_GT)
