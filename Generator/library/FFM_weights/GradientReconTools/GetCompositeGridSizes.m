function [n_comp, n_comp_pad] = GetCompositeGridSizes(n_ap, n_WFS_frames, wind_vecs, ssp)
%
%  n_comp = GetCompositeGridSizes(n_ap, n_WFS_frames, wind_vecs);
%  
% This function computes the grid sizes needed to store the composite
% high resolution gradients and phases for a set of FFM frames.
%
% Input:
%   n_ap         - number of pixels across the aperture (diameter)
%   n_WFS_frames - vector containing the number of frames used for each
%                  FFM reconstruction
%   wind_vecs    - 2D double array that specifies angle (theta) and 
%                  magnitude (r) of wind velocity for each layer:
%                      wind_vecs(L,1) = r   
%                      wind_vecs(L,2) = theta 
%                  L = 1, 2, ..., n_layers.
%       ssp      - subsample parameter
%                          Modified November, 2015, to allow non square 
%                          subapertures. (J. Nagy) So now, ssp can
%                          be a scalar (equal diameters in both directions)
%                          or a 1-by-2 array, [ssp_rows, ssp_cols]
%
% Output:
%   n_comp     - vector containing composite grid sizes
%   n_comp_pad - amount the aperture needs to be padded to get to the
%                composite grid sizes
%

%
%  J. Nagy
%  October, 2013
%

%
%  Convert wind_vecs into pixel shif\t values:
%
%
%  Modified November, 2015, to allow non square subapertures.
%  J. Nagy
%
%
%  Modified November, 2015, to allow non square subapertures.
%  J. Nagy
%
if length(ssp) == 1
    ssp = [ssp, ssp];
end
deltax = ssp(2)*wind_vecs(:,1) .* cos(wind_vecs(:,2));
deltay = ssp(1)*wind_vecs(:,1) .* sin(wind_vecs(:,2));
nframes = length(n_WFS_frames);
n_comp = zeros(nframes,1);
for k = 1:nframes
  %
  %  We need to find a composite grid size, and the amount of padding needed
  %  to get from n_ap to n_conmp.  
  %
  [n_comp(k), n_comp_pad] = GetCompositeGridSize(n_ap, n_WFS_frames(k), deltax, deltay);
end