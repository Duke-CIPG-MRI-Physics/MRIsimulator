function [regionA] = getRegionA(Matrix_Size_Acquired,FOV_acquired,pA,PF_Factor,R)
%Calculates region A for TWIST
%If pA = .05, n_pixels_in_A / n_pixels_acquired = .05

%% CHEATING CODE
%Please remove this section once GRAPPA and PF have been properly
%implemented

PF_Factor = [6/8 6/8];
R = [2 3];


%%
%first we must calculate n_pixels_acquired based on PF
if length(PF_Factor)>1
    Undersampled_Matrix_Size = [round(Matrix_Size_Acquired(2) * PF_Factor(1)), round(Matrix_Size_Acquired(3) * PF_Factor(2))];
    n_pixels_acquired = Undersampled_Matrix_Size(1)*Undersampled_Matrix_Size(2);
else
    n_pixels_acquired = Matrix_Size_Acquired(2)*Matrix_Size_Acquired(3);
end

%next we account for GRAPPA
if length(R)>1
    n_pixels_acquired = n_pixels_acquired/(prod(R));
end

%next we calculate how many pixels we will be included in A
n_pixels_in_A = round(pA * n_pixels_acquired);

%we want region A to be isotropic in spatial resolution, we must define a
%grid of k-space frequenices to do this

kspace_pixel_size = 1./FOV_acquired; 
kspace_extent = kspace_pixel_size .* (Matrix_Size_Acquired-1)/2;

phase_frequencies = -kspace_extent(2):kspace_pixel_size(2):kspace_extent(2);
slice_frequencies = (-kspace_extent(3):kspace_pixel_size(3):kspace_extent(3))';

[phase_mesh,slice_mesh] = meshgrid(slice_frequencies,phase_frequencies);
frequency_grid = sqrt(phase_mesh.^2+slice_mesh.^2);
frequency_grid_idx = reshape(1:numel(frequency_grid), size(frequency_grid));

%We create a table which tracks each pixel's frequency and index
frequency_table = table(frequency_grid_idx(:),frequency_grid(:),'VariableNames',{'Linear Index','Frequency'});
%We sort the table from low to high frequencies
frequency_table = sortrows(frequency_table,2,'ascend');
%We keep the first n_pixels_in_A number of pixels from the sorted list
regionA_table = frequency_table{1:n_pixels_in_A,:};
%Those points now define region A
regionA = zeros(size(frequency_grid));
regionA(regionA_table(:,1)) = 1;
regionA = logical(regionA);

figure
imshow(regionA)
title('Region A within Acquired Matrix')

end