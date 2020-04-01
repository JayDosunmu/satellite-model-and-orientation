function image_frames = MosaicPadding(image_frames, pupil_mask)
%
%  When building the mosaic images from set of frames, there
%  might be issues from boundary artifacts if the region of interest
%  is not known throughout the whole domain.  In particular, it
%  it is only known in a region defined by the pupil_mask.
%  This function will use inpainting techniques to fill in the
%  missing regions, and hopefully avoid the boundary artifacts
%  when we use ReconMosaic.m
%
%  Note that if we are using this for a multi-layer FFH problem,
%  we might build mosaics for each layer.  But we need only do
%  this padding once.
%

[m, n, nframes] = size(image_frames);

%h = waitbar(0,'Start padding ...');
for k = 1:nframes
  X = image_frames(:,:,k);
  X(~pupil_mask) = NaN;
  image_frames(:,:,k) = inpaintn(X);
  %waitbar(k/nframes,h)
end
%close(h)
