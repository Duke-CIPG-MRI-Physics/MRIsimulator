%%
clear
close all
clc
%% Load file
filename = uigetfile('*.*','Select the Segmentation File');

underscores = find(filename == '_');

Dimension = [str2double(filename(underscores(2)+1:underscores(3)-1)); ...
            str2double(filename(underscores(3)+1:underscores(4)-1));...
            str2double(filename(underscores(4)+1:underscores(5)-1))];
        
Type = filename(underscores(5)+1:underscores(6)-1);

BreastData = uint8(ReadRaw(filename,Dimension(1),Dimension(2),...
    Dimension(3),Type,'LittleEndian'));
%% Change orientation
NewOrientation = permute(BreastData,[3,1,2]); % axial view

% NewOrientation = permute(BreastData,[3,2,1]); % sagittal view

NewDimension = size(NewOrientation);
%% Save as .raw file with new name reflective of new orientation
WriteRaw(strcat(filename(1:underscores(2)),num2str(NewDimension(1)),...
    '_',num2str(NewDimension(2)),'_',num2str(NewDimension(3)),'_',Type,...
    '_LE.raw'),NewOrientation,Type,'LittleEndian')
%% Save as .hdr/.img with new name reflective of new orientation
DualFileStructure = make_nii(NewOrientation);

save_nii(DualFileStructure,strcat(filename(1:underscores(2)),...
    num2str(NewDimension(1)),'_',num2str(NewDimension(2)),...
    '_',num2str(NewDimension(3)),'_',Type,'_LE'));