function imageWithPhase = addPhaseToImage(inputImage)
% Get image size
[Nx, Ny, Nz] = size(inputImage);

%Generate smooth quadratic phase
[x,y,z] = ndgrid(linspace(-0.5,0.5,Nx), ...
                 linspace(-0.5,0.5,Ny), ...
                 linspace(-0.5,0.5,Nz));
phi     = 2*pi*( 0.4*(x.^2 + y.^2) + 0.2*z.^2 );     % radians
imageWithPhase = inputImage .* exp(1i*phi);
end