function psfs_poly  = p_dora_color_psf(phase, PSF_params)
%
%  This function constructs a polychromatic PSF, which is an average of
%  monochromatic PSFs at a given set of wavelengths.  Each of the 
%  monchformatic PSFs uses the standard Fourier optics model:
%
%        PSF(lambda) = ifft2(pupil*exp(i*phase))
%
%  05/29/2017 - The phase argument are the phases in a FFM

color_planes = PSF_params.color_planes;
wavelength   = PSF_params.wavelength;
bandwidth    = PSF_params.bandwidth;
PupilDiam    = PSF_params.PupilDiam;
n_obj        = PSF_params.n_obj;
pupil_poly   = PSF_params.pupil_amp;

pupild = p_GetPupilDiams(wavelength, bandwidth, PupilDiam, color_planes);

h = zeros(n_obj, n_obj, size(phase,3), length(pupild));
for k = 1:length(pupild)
  poly_scale = PupilDiam/pupild(k);
  for j = 1:size(phase,3)
    phi  = phase(:,:,j)./poly_scale;
    P    = pupil_poly.*exp(sqrt(-1)*phi);
    psf1 = ifftshift(abs(ifft2(P)).^2);
    %
    % The next function does an interpolation.  Since we always use
    % 'same' in the original code, I've changed this to eliminate the
    % input variable 'same'.
    %
    % psf1 = ImgScale(psf1, poly_scale, 'same');
    psf1 = p_ImgScale(psf1, poly_scale);
    h(:,:,j,k) = psf1;
  end        
end
psfs_poly = squeeze(mean(h,3));
