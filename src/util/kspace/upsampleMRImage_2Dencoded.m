function [upsampledIm] = upsampleMRImage_2Dencoded(im, upsampled2DSize)
% This function just zeropads k-space to "upsample" an image while
% preserving the measured k-space data. This function correctly handles
% odd and even matrix sizes of both the input image and upsampled Size. the
% input image can be n-dimmensional, but the first two dimmensions will be
% interpreted as the 2D encoded dimmension. The input upsampledSize must be
% a vector of size [m n] where m and n are the size of the output image.
% The output upsampledIm will have size [m n imSize(3:end)]
%
% Usage:
% upsampledIm = upsampleMRImage_2Dencoded(im, upsampledSize)

inputSz = size(im);
input2DSize = inputSz(1:2);
outputSz = inputSz;
outputSz(1:2) = upsampled2DSize;

% Changing the size will also add a scaling factor. We will remove it to
% preserve image intensities.
intensityScaleFactor = prod(upsampled2DSize./input2DSize);

% Take 2D FFT of each slice to get k-space image with the DC frequency in
% the center of k-space
centeredKsp = fftshift(fftshift(fft2(im),1),2);

% Calculate location of old k-space center pixel
originalKSpaceCenter = floor(0.5*input2DSize) + 1;

% Calculate the new center pixel of k-space
paddedKSpaceCenter = floor(0.5*upsampled2DSize) + 1;

% Zeropad upsampled k-space
paddedKspace = complex(zeros(outputSz));
startVal = paddedKSpaceCenter-originalKSpaceCenter + 1;
endVal = paddedKSpaceCenter-originalKSpaceCenter + input2DSize;
paddedKspace(startVal(1):endVal(1),startVal(2):endVal(2),:,:,:,:,:) = centeredKsp;

% Take 2D IFFT of each slice to get upsampled image
upsampledIm = ifft2(ifftshift(ifftshift(paddedKspace,1),2));

%Undo intensity scaling
upsampledIm = upsampledIm*intensityScaleFactor;

end