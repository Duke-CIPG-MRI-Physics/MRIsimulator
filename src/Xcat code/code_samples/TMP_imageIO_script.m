

imageFile = 'E:\BreastSegmentationNew\CTA0296_SEGMENT_451_470_297_uint8_LE.raw';
Raw = ReadRaw(imageFile,451,470,297,'int8');



% MatLab's internal orientation is different that ImageJ if you want the
% same orientation when you use imshow or imtool you can use permute and
% flipdim just be sure to undo this when outputting back out.


% Set internal 6's to 5's.
imageFile = 'E:\BreastSegmentationNew\CTA0296_SKINMASK_451_470_297_uint8_LE.raw';
Skin = ReadRaw(imageFile,451,470,297,'int8');

Raw(Raw>5 & Skin==0) = 5;

WriteRaw('E:\BreastSegmentationNew\CTA0296_SEGMENT_451_470_297_uint8_LE_PROCESSED.raw',Raw,'int8')

