function pupild = p_GetPupilDiams(lambda0, dlambda, pupil0, nColors)
%
%  When building a PSF, we need the pupil aperture function.  In the case 
%  of a polychromatic PSF, the pupil aperture size (i.e., its diameter)
%  depends on the wavelength.  This function gets the pupil diameters for
%  the different wavelengths.
%
%  Input:  lambda0 = starting wavelength of observation in microns 
%                    (wavelength of phvars)
%                    For example, 0.750
%          dlambda = optical bandwidth
%                    For example, 0.3
%          pupil0  = diameter of pupil at wavelength lambda0
%                    For example, 192
%          nColors = number of wavelengths.
%                    For example, 13
%
%  Output:  pupild = vector of length nColors containing the pupil
%                    diameters.
%

if dlambda>0
    pupil_max = (1 + dlambda/(2*lambda0))*pupil0;
    pupil_max = 2*round(pupil_max/2);

    pupil_min = (1 - dlambda/(2*lambda0))*pupil0;
    pupil_min = 2*round(pupil_min/2);

    %pupild = pupil_min:2:pupil_max;
    pupild = pupil_min:pupil_max;

    idx0 = find(pupild == pupil0);

    idx = round(linspace(1, length(pupild), nColors));
    if ~any(idx(:) == idx0)
      %  In this case, the pupil diameter spacing strange -- trying to fix
      %  I'm not sure if this conditional statement is needede.  Might be
      %  good to remove.
    %  idx(round(length(idx)/2)) = idx0;
    end
    pupild = pupild(idx);
else     
    pupild = pupil0;
end



