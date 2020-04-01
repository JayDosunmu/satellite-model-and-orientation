function pupil_maskC = GetCompositeMask(pupil_mask,wind_vecs, ssp, n_ap, n_comp, n_comp_pad, n_frames, n_layers)
%
%  pupil_maskC = GetCompositeMask(pupil_mask, wind_vecs, ssp, n_ap, ...
%        n_comp, n_comp_pad, n_frames, n_layers);
%
%  This function comptues the mask of the composit grid.
%
%  Input:
%    pupil_mask  - pupil mask
%    wind_vecs   - specifies angle (theta) and magnitude (r) of wind 
%                  velocity for each layer:
%                       wind_vecs(L,1) = r   
%                       wind_vecs(L,2) = theta 
%                  L = 1, 2, ..., n_layers.
%    ssp         - subsampling parameter that relates low resolution
%                  grid to high resolution grid
%    n_ap        - size of aperture
%    n_comp      - size of composite grid
%    n_comp_pad  - size of padding to get from n_ap to n_comp
%                  (see GetCompositeGridSize.m)
%    n_frames    - number of frames of data
%    n_layers    - number of atmospheric layers
%
%  Output:
%    pupil_maskC - 3D logical array containing mask for composite
%                  information (phases and gradients) on leach layer.
%

%
%  J. Nagy
%  August, 2012
%

%
%    * Start with the pupil mask on the composite grid.
%    * Move the pupil mask across the composite grid
%    * BUT make sure that the pupil mask moves in the opposite 
%      direction as the sky.  This done by adding pi to the
%      angle of the wind vectors.
%
deltax2 = ssp*wind_vecs(:,1) .* cos(pi+wind_vecs(:,2));
deltay2 = ssp*wind_vecs(:,1) .* sin(pi+wind_vecs(:,2));

A2 = MotionMatrix(deltax2, deltay2, n_comp, n_ap, n_frames);
pupil_mask_pad = padarray(pupil_mask, [n_comp_pad,n_comp_pad],0,'both');
N = n_comp*n_comp;
v = kron(speye(n_layers),pupil_mask_pad(:));
v = full(A2*v);
v = reshape(v, N, n_frames, n_layers);
pupil_maskC = reshape(squeeze(sum(v, 2) > 0),n_comp,n_comp,n_layers);
