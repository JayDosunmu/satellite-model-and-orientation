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
deltax = ssp*wind_vecs(:,1) .* cos(wind_vecs(:,2));
deltay = ssp*wind_vecs(:,1) .* sin(wind_vecs(:,2));
nframes = length(n_WFS_frames);
n_comp = zeros(nframes,1);
for k = 1:nframes
  %
  %  We need to find a composite grid size, and the amount of padding needed
  %  to get from n_ap to n_conmp.  
  %
  [n_comp(k), n_comp_pad] = GetCompositeGridSize(n_ap, n_WFS_frames(k), deltax, deltay);
end