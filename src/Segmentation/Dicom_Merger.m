%Roberto Carrasscosa 1/14/25
% 
%Designed for converting non dicom volume to dicom series with given
%header info
%
%Takes input volume and merges with header info to create dicom folder
%
%Accepts a .dcm or structure file for header file and reads header info
%
%Accepts .nii and .tiff files for volume
%
%Creates a new folder in working directory with same name as volume
%file and saves files there

clear;clc;
Header_File = input("Header File Name: ");
Volume_File = input("Volume File Name: ");

%Destination Folder Setup
folder_name = strcat("Fixed_",Volume_File,".dcm");

if exist(folder_name, 'dir') % Check if the folder already exists
    error("Desintation Folder Already Exists!")
end

%Reading header file
if isstruct(Header_File) == 1 %Check if structure file
    header = Header_file;
elseif endsWith(Header_File,".dcm")
    header = dicominfo(Header_File); %read dicominfo if not
else
    error("Header Filetype not Supported!")
end

%Reading volume file
if endsWith(Volume_File,".nii") %Check for niftii
    volume = niftiread(Volume_File); %rotate to correct orientation
    volume = rot90(volume,3);
    volume = fliplr(volume);
elseif endsWith(Volume_File,".tiff") || endsWith(Volume_File,".tif") %check for tiff
    volume = tiffreadVolume(Volume_File);
else
    error("Volume Filetype not Supported!")
end

%Ask if user wants preview prior to writing

preview_choice = input("Preview Volume File prior to writing? (y/n): ",'s');
z_height = size(volume,3);

if preview_choice == "y"
    imshow(volume(:,:,ceil(z_height/2)),[])
    preview_acceptance = input("Continue to write? (y/n): ",'s');

    if preview_acceptance == "y"
        write_confirm = 1;
    elseif preview_acceptance =="n"
        write_confirm = 0;
    end

elseif preview_choice == "n"
    write_confirm = 1;
end

if write_confirm == 1
mkdir(folder_name); % Create the folder
    %Merging
for ii = 1:z_height
    filename = sprintf("Fibroglandular_Mask_V1_%d",ii);
    file_path = fullfile(folder_name,filename);
    dicomwrite(volume(:,:,ii),file_path,header); %Write volume and header to dicom series in folder we created
end    
end



