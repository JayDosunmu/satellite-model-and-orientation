function phase = p_ReconstructPhase(phase_gradx, phase_grady, pupil_mask, scale, alpha)
%
%    phase = p_ReconstructPhase(phase_gradx, phase_grady, pupil_mask, alpha);
%
%  This function reconstructs the phase from its gradient measurements.
%  We assume a Fried geometry of the wavefront sensor.
%
%  The p is used to denote that this MATLAB code has been stripped of 
%  any unnecessary computations to try to get it ready for the parallel
%  implementation.
%
%  Input:
%    phase_gradx - array containing x-gradients of phases.  This can be
%                  either a 2D or a 3D array.  In case of 3D, each 2D
%                  slice, phase_gradx(:,:,k), should be gradients for a 
%                  particular frame of data.
%    phase_grady - array containing y-gradients of phases.  This can be
%                  either a 2D or a 3D array.  In case of 3D, each 2D
%                  slice, phase_gradx(:,:,k), should be gradients for a 
%                  particular frame of data.
%    pupil_mask  - pupil mask for phases.  Must be a 2D array, so
%                  that each phase has the same pupil mask.
%
%  Optional Input:
%    scale       - allows for a normalization scale for specific geometries.
%                  Default is to use scale = 1;
%    alpha       - Tikhonov regularization parameter. Default is
%                  alpha = 1e-6.
%
%  Output:
%    phase       - reconstructed phases.
%

%
%  J. Nagy
%  November, 2013
%

%
%  Check inputs and set default values.
%
switch nargin
  case 3
    scale = []; alpha = []; 
  case 4
    alpha = [];
end

if isempty(scale)
  scale = 1;
end
if isempty(alpha)
  %
  %  Note that the default for this should be different in the C-code
  %  because we are using a different method to solve the regularized
  %  system.
  %
  alpha = 1e-6;
end

alpha = scale*alpha;

%
%  This can help reduce some boundary artifacts ...
%
pad_size = 1;
phase_gradx = padarray(phase_gradx, [pad_size, pad_size], 'both');
phase_grady = padarray(phase_grady, [pad_size, pad_size], 'both');
pupil_mask  = padarray(pupil_mask, [pad_size, pad_size], 'both');

%
%  Get dimensions of padded versions ...
%
[m_ap, n_ap, n_frames] = size(phase_gradx);

%
%  Now set pu matrices and solve ...
%
Dx = GradxMatrix(n_ap, scale);
Dy = GradyMatrix(n_ap, scale);

P = spdiags(pupil_mask(:), 0, n_ap*n_ap, n_ap*n_ap);
A = [P*Dx; P*Dy; alpha*speye(n_ap*n_ap)];
b = [reshape(phase_gradx,n_ap*n_ap,n_frames); reshape(phase_grady,n_ap*n_ap, n_frames); zeros(n_ap*n_ap,n_frames)];
phase_vecs = A\b;
pupil_mask = pupil_mask(:,:,ones(1,n_frames));
phase = pupil_mask .* reshape(phase_vecs, n_ap, n_ap, n_frames);

%
%  Extract phase from padded version ...
%
phase = phase(pad_size+1:m_ap-pad_size, pad_size+1:n_ap-pad_size, :);
