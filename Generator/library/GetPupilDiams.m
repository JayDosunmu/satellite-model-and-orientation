function pupild = GetPupilDiams(lambda0, dlambda, pupil0, nColors)
%
%  lambda0 = starting wavelength of observation in microns (wavelength of phvars)
%            Default: 0.750
%  dlambda = optical bandwidth
%            Default: 0.3
%  pupil0  = diameter of pupil at wavelength lambda0
%            Default: 192
%  nColors = number of wavelengths.
%            Default: 13
%

%
% Check input parameters, and set default values if necessary:
%
if nargin == 0, lambda0=[]; dlambda=[]; pupil0=[]; nColors=[]; end
if nargin == 1, dlambda=[]; pupil0=[]; nColors=[]; end
if nargin == 2, pupil0=[]; nColors=[]; end
if nargin == 3, nColors=[]; end
if isempty(lambda0), lambda0 = 0.75; end
if isempty(dlambda), dlambda = 0.3; end
if isempty(pupil0), pupil0 = 192; end
if isempty(nColors), nColors = 13; end

pupil_max = (1+dlambda/(2*lambda0))*pupil0;


if mod(pupil0,2)==0
    pupil_max = 2*round(pupil_max/2);
else
    pupil_max = 2*round(pupil_max/2)+1;
end

pupil_min = (1 - dlambda/(2*lambda0))*pupil0;
if mod(pupil0,2)==0
    pupil_min = 2*round(pupil_min/2);
else
    pupil_min = 2*round(pupil_min/2)+1;
end
% 
% 
% pupil_max = 2*round(pupil_max/2);
% 
% pupil_min = (1-dlambda/(2*lambda0))*pupil0;
% pupil_min = 2*round(pupil_min/2);

pupild = pupil_min:2:pupil_max;

idx0 = find(pupild == pupil0);
% pupild = 125;
return

idx = round(linspace(1, length(pupild), nColors));
if ~any(idx(:) == idx0)
%  warning('pupil diameter spacing strange -- trying to fix')
  idx(round(length(idx)/2)) = idx0;
end
pupild = pupild(idx);




