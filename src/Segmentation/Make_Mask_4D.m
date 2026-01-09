%Roberto Carrascosa

%Used for converting several 3D .tif mask files into a single 4D logical volume
%file for ease of use later

clear;clc;

Vol_1 = logical(tiffreadVolume("Fat Mask_V1.tif"));
Vol_2 = logical(tiffreadVolume("Fibroglandular Mask_V1.tif"));
Vol_3 = logical(tiffreadVolume("Lesion Mask_V1.tif"));

List_of_Masks = ["Fat" "Fibroglandular" "Lesion"];

Combined_Masks_4D = cat(4,Vol_1,Vol_2,Vol_3);

%save('V1_Maks_4D','Combined_Masks_4D','List_of_Masks')