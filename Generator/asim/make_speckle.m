%
%  Code to generate a speckle image with a given D/r0
%
%  Oversampling increases the sampling rate by this factor
%  above Nyquist, and the makes speckles are easier to see,
%  although it is not efficient for processing.
%
% function [speckle] = make_speckle(D, r0, oversampling, dim_aperture)
%
%  Inputs       dim_aperture    - the number of sample points for the aperture
%           D - Telescope diameter corresponding to the number of points
%           r0 - Fried's parameter (vector) for computing PSF's integrated
%           over spectral band
%           oversampling - the increase over Nyquist of the speckle
%                                                                       sampling
%
%
%  Outputs      speckle - a speckle image of size dim_aperture*oversampling
%                                                               normalised to sum to 1.
%
%  Uses genphase.m, kolmogorov.m
%
%   Richard Lane June 2001
%   Department of Electrical and Electronic Engineering
%   University of Canterbury
%   Christchurch
%   New Zealand
%
%   Code may be copied and used freely provided acknowledgement of the original source
%   is maintained.
%
function [int_speckle,phase] = make_speckle(D, r0, oversampling, dim_aperture, photons)
%
aperture = circle2(dim_aperture/2, dim_aperture);

nB = length(r0);
int_speckle = zeros(2*dim_aperture*oversampling);


phase2 = kolmogorov(dim_aperture, D./r0);

for b=1:nB
    phase       = phase2(:,:,b);
    cmplx_aper  = exp(sqrt(-1)*phase).*aperture;
    cmplx_aperture = zeros(2*dim_aperture*oversampling);
    cmplx_aperture(1:dim_aperture,1:dim_aperture) = cmplx_aper;
    cmplx_image = ifft2(cmplx_aperture)*(2*dim_aperture);
    speckle = fftshift(abs(cmplx_image).^2);

    speckle = speckle/sum(sum(speckle));
    int_speckle = int_speckle + speckle;
end
