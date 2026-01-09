function [coilImgs, coilKspace, coilSens] = ...
    simulateCoils(inputVol, nCoilsEven, sigma, coilSplitDirection, coilDirection, opts)
%SIMULATECOILS  Synthetic multi-coil simulation for 3-D MRI data.
%
%   [coilImgs, coilKspace, coilSens, coilCenters] = ...
%       simulateCoils(inputVol, nCoilsEven, sigma, ...
%                     coilSplitDirection, coilDirection, opts, bTest)
%
%   INPUT ARGUMENTS
%   ------------------------------------------------------------------
%   inputVol           :  Nx-by-Ny-by-Nz real or complex array
%
%   nCoilsEven         :  even integer (>=2). Elements mirrored on two
%                         planes ±planeOffset along coilSplitDirection.
%
%   sigma              :  Gaussian widths (normalised FOV units, -1…+1)
%                         scalar      – same σ in x,y,z
%                         [σx σyz]    – anisotropic (σy=σz)
%                         [σx σy σz]  – fully anisotropic
%
%   coilSplitDirection :  'x' | 'y' | 'z'   (axis of the two half-arrays)
%   coilDirection      :  'x' | 'y' | 'z'   (orthogonal axis to space coils)
%
%   opts (struct, *all fields optional*)
%     .planeOffset  = 1.10    distance of coil planes (>|1|)
%     .phaseAmp     = 1.0     amplitude of linear phase (radians)
%     .rngSeed      = 0       0⇒deterministic; else seed for RNG
%     .dtype        = "single"|"double"   storage precision  (default "single")
%     .noiseRMS     = 0       complex RMS of k-space noise per point, per coil
%                             (set >0 to add noise)
%
%   OUTPUTS
%   ------------------------------------------------------------------
%   coilImgs     :  Nx×Ny×Nz×nCoils  complex  (image-space per coil)
%   coilKspace   :  Nx×Ny×Nz×nCoils  complex  (k-space per coil, noise added)
%   coilSens     :  Nx×Ny×Nz×nCoils  complex  (normalised sensitivities)
%
%   ------------------------------------------------------------------
%   Scott H. Robertson  —  July 2025
%   ------------------------------------------------------------------

% ------------------ 0. INPUT VALIDATION ---------------------------------
arguments
    inputVol {mustBeNumeric, mustBeNonempty}
    nCoilsEven (1,1) {mustBeInteger, mustBePositive}
    sigma {mustBeNumeric, mustBeVector, mustBeNonempty}

    coilSplitDirection (1,1) char  {mustBeMember(coilSplitDirection,{'x','y','z'})}
    coilDirection      (1,1) char  {mustBeMember(coilDirection,     {'x','y','z'})}

    opts.planeOffset   (1,1) double  {mustBePositive}    = 1.10
    opts.phaseAmp      (1,1) double  {mustBeNonnegative} = 1.0
    opts.rngSeed       (1,1) double                      = 0
    opts.dtype         (1,1) string {mustBeMember(opts.dtype,["single","double"])} ...
                                                        = "single"
    opts.noiseRMS      (1,1) double  {mustBeNonnegative} = 100
end

if mod(nCoilsEven,2)~=0
    error("nCoilsEven must be even (mirrored halves).");
end
if coilSplitDirection==coilDirection
    error("coilSplitDirection and coilDirection must differ.");
end

% σ interpretation
sigma = double(sigma(:));
switch numel(sigma)
    case 1, sigma = [sigma sigma sigma];
    case 2, sigma = [sigma(1) sigma(2) sigma(2)];
    case 3  % ok
    otherwise
        error("sigma must be scalar, 2-element, or 3-element vector.");
end
[sx,sy,sz] = deal(sigma(1),sigma(2),sigma(3));

% RNG seed
if opts.rngSeed~=0, rng(opts.rngSeed,'twister'); end

% axis index helper
ax = struct('x',1,'y',2,'z',3);
splitAx  = ax.(coilSplitDirection);
spreadAx = ax.(coilDirection);

% ------------------ 1. COORDINATE GRID ----------------------------------
[Nx,Ny,Nz] = size(inputVol);
dtypeFun   = str2func(char(opts.dtype));

[x,y,z] = ndgrid( linspace(-1,1,Nx), ...
                  linspace(-1,1,Ny), ...
                  linspace(-1,1,Nz) );
x=dtypeFun(x); y=dtypeFun(y); z=dtypeFun(z);

% ------------------ 2. COIL CENTRES -------------------------------------
nSide = nCoilsEven/2;
spread = linspace(-1,1,nSide).';          % nSide×1

c1 = zeros(nSide,3);  c2 = c1;            % negative / positive planes
c1(:,splitAx) = -opts.planeOffset;
c2(:,splitAx) =  opts.planeOffset;
c1(:,spreadAx)= spread;
c2(:,spreadAx)= spread;

coilCenters = [c1;c2];                    % nCoilsEven×3

% ------------------ 3. PRE-ALLOCATE -------------------------------------
coilSens   = zeros(Nx,Ny,Nz,nCoilsEven, char(opts.dtype));
coilImgs   = zeros(size(coilSens),        char(opts.dtype));
coilKspace = zeros(size(coilSens),        char(opts.dtype));

% ------------------ 4. SENSITIVITIES ------------------------------------
for c = 1:nCoilsEven
    xc = coilCenters(c,1); yc = coilCenters(c,2); zc = coilCenters(c,3);

    % magnitude (Gaussian)
    Smag = exp( -0.5*((x-xc)/sx).^2 ...
                -0.5*((y-yc)/sy).^2 ...
                -0.5*((z-zc)/sz).^2 );

    % phase (linear)
    theta = 2*pi*(c-1)/nCoilsEven;
    a = opts.phaseAmp*cos(theta);
    b = opts.phaseAmp*sin(theta);
    cph = 0.3*opts.phaseAmp*sin(theta);    % small z term
    phi = a*x + b*y + cph*z;

    coilSens(:,:,:,c) = dtypeFun( Smag .* exp(1i*phi) );
end

% normalise so Σ|S|² = 1
normFactor = sqrt(sum(abs(coilSens).^2,4) + eps(dtypeFun(1)));
coilSens   = coilSens ./ normFactor;

% ------------------ 5. COIL IMAGES & k-SPACE ----------------------------
for c = 1:nCoilsEven
    coilImgs(:,:,:,c)   = inputVol .* coilSens(:,:,:,c);
    coilKspace(:,:,:,c) = fftshift( fftn(coilImgs(:,:,:,c) ) );
end

% ------------------ 6. ADD COMPLEX WHITE NOISE (optional) ---------------
if opts.noiseRMS > 0
    nRMS = opts.noiseRMS / sqrt(2);      % std-dev per real/imag part

    % choose float precision consistent with coilKspace
    noiseClass = class(coilKspace);

    for c = 1:nCoilsEven
        % generate noise just for this coil (Nx×Ny×Nz)
        noise = nRMS * ( randn(size(coilKspace(:,:,:,c)), noiseClass) + ...
                         1i * randn(size(coilKspace(:,:,:,c)), noiseClass) );

        % add in-place
        coilKspace(:,:,:,c) = coilKspace(:,:,:,c) + noise;

        % regenerate the noisy image *only* for this coil
        coilImgs(:,:,:,c)   = ifftn( ifftshift(coilKspace(:,:,:,c) ) );
    end
end
end

