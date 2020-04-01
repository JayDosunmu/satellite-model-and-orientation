function pupil_maskC = GetCompositeMask(A, pupil_mask, n_comp, n_comp_pad, n_frames, n_layers)
%
%  pupil_maskC = GetCompositeMask(pupil_mask, wind_vecs, ssp, n_ap, ...
%        n_comp, n_comp_pad, n_frames, n_layers);
%
%  This function comptues the mask of the composit grid.
%
%  Input:
%    A           - matrix that defines the motion of the various layers
%    pupil_mask  - pupil mask
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
%  Modified April, 2016 to reduce memory requirements.
%

%
%    * Start with the pupil mask on the composite grid.
%    * Move the pupil mask across the composite grid
%    * BUT make sure that the pupil mask moves in the opposite 
%      direction as the sky.  This done by adding pi to the
%      angle of the wind vectors.
%
%if length(ssp) == 1
%    ssp = [ssp, ssp];
%end
%deltax2 = ssp(2)*wind_vecs(:,1) .* cos(pi+wind_vecs(:,2));
%deltay2 = ssp(1)*wind_vecs(:,1) .* sin(pi+wind_vecs(:,2));
%A2 = MotionMatrix(deltax2, deltay2, n_comp, n_ap, n_frames);
pupil_mask_pad = padarray(pupil_mask, [n_comp_pad,n_comp_pad],'both');
N = n_comp*n_comp;
v = kron(speye(n_layers),pupil_mask_pad(:));
%v = full(A2*v);
%
%  Instead of creating a new matrix that does the motion in the opposite
%  direction, let's move in the given direction, then rotate the result.
%  This will avoid having to create another large matrix.
%
v = full(A*v);
v = reshape(v, N, n_frames, n_layers);
for k = 1:n_layers
%   v(:,:,k) = transpose(rot90(transpose(v),2));
   v(:,:,k) = transpose(rot90(transpose(v(:,:,k)),2));

end
pupil_maskC = reshape(squeeze(sum(v, 2) > 0),n_comp,n_comp,n_layers);
