function [intensities,timepoints] = Generate_Intensity_Plots(output_phantom,temporal_resolution,position)
% Generates plot of intensity at a given position in phantom over time


arguments (Input)
    output_phantom (:,:,:,:)
    temporal_resolution (1,1) 
    position (1,3)
end

intensities = squeeze(output_phantom(position(1),position(2),position(3),:));
timepoints = prep_scan_time:temporal_resolution:temporal_resolution*size(output_phantom,4);


end