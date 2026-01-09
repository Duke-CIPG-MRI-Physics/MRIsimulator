function [twistMasks,Tframe_sec,Trecon_sec] = calculateTwistMasks(FOV_m,matrixSize,A,B,TR_sec,R,partialFourier)
% Calculate the sampling mask considering the parallel imaging factors
% and partial fourier
tempCalRegion = 24*ones(size(R)); % The calibration region size doesnt matter for these simulations
[kspace_mask, ~] = generate_sampling_masks(matrixSize, R, partialFourier, tempCalRegion,false);

% Calculate the number of phase encodes in A and B regions
nPE_total = sum(kspace_mask(:)); % Total number of phase encodes
nPE_Aframe = ceil(nPE_total.*A);

% Calculate the k-space distance from DC frequency, accounting for non isotropic voxels
[kx, ky, kz] = calculate_spatial_frequencies(FOV_m, matrixSize);
[Mz,My] = meshgrid(kz,ky);
Mr = sqrt(My.^2 + Mz.^2);
Mr(~kspace_mask) = inf;

% There may be a little shifting of the number of A phase encodes...
[sortR, sortI] = sort(Mr(:),'ascend');
MaxR_A = sortR(nPE_Aframe);
mask_A = (Mr <= MaxR_A);
nPE_Aframe = sum(mask_A(:));

% Calculate the number of B phase encodes in each recon and frame
nPE_Brecon = nPE_total - nPE_Aframe;
nPE_Bframe = ceil(nPE_Brecon.*B);

% Calculate the total number of frames to reconstruct one TWIST time
% point
nFramesPerRecon = ceil(nPE_Brecon/nPE_Bframe);

% Calculate B mask
mask_B = double(kspace_mask & (~mask_A));
mask_B_idx = find(mask_B);
mask_B_idx = mask_B_idx(randperm(length(mask_B_idx))) % randomize indices
for iB = 1:nFramesPerRecon
    % take first nPE_Bframe of randomized Idx
    rndIdxLookup = (iB-1)*nPE_Bframe+(1:nPE_Bframe);
    if(any(rndIdxLookup > length(mask_B_idx)))
        rndIdxLookup(rndIdxLookup > length(mask_B_idx)) = [];
    end
    frameIdx = mask_B_idx(rndIdxLookup);

    % Save the frame number (starting at 2 since 1 is the A region)
    mask_B(frameIdx) = 1+iB;
end

% Calculate the number of phase encodes per frame
nPE_singleFrame = nPE_Aframe + nPE_Bframe;

% Calculate the time required to acquire data for a single frame
Tframe_sec = nPE_singleFrame * TR_sec;

% Calculate the time required to reconstruct one Twist time point
Trecon_sec = Tframe_sec*nFramesPerRecon;

twistMasks = mask_A + mask_B;
end