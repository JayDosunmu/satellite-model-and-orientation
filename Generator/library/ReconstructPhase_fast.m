function phase = ReconstructPhase_fast(phase_gradx, phase_grady, pupil_mask, scale, alpha, method)
%
%    phase = ReconstructPhase_fast(phase_gradx, phase_grady, pupil_mask, alpha);
% 
%  This function reconstructs the phase from its gradient measurements.
%  We assume a Fried geometry of the wavefront sensor.
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
%    method      - character string, either 'direct' or 'lsqr'.
%                  direct:    Tikhonov regularized least squares problems
%                             are solved by MATLAB's backslash operator,  
%                             exploiting sparse matrix tools.
%                  iterative: Tikhonov regularized least squares problems
%                             are solved using the iterative method LSQR.
%                  Default is 'direct'.
%
%  Output:
%    phase       - reconstructed phases.
%

%
%  J. Nagy
%  November, 2012
%

%
%  Check inputs and set default values.
%
% narginchk(3, 6);
switch nargin
  case 3
    scale = []; alpha = []; method = [];
  case 4
    alpha = []; method = [];
  case 5
    method = [];
end

if isempty(scale)
  scale = 1;
end
if isempty(alpha)
  alpha = 1e-6;
end
if isempty(method)
  method = 'direct';
end
alpha = scale*alpha;

%
%  This can help reduce some boundary artifacts ...
%
pad_size = 1;
phase_gradx = padarray(phase_gradx, [pad_size, pad_size],0, 'both');
phase_grady = padarray(phase_grady, [pad_size, pad_size],0, 'both');
pupil_mask = padarray(pupil_mask, [pad_size, pad_size],0, 'both');

[m_ap, n_ap, n_frames] = size(phase_gradx);

Dx = GradxMatrix(n_ap, scale);
Dy = GradyMatrix(n_ap, scale);

% For this code, let's assume the pupil maskes are the same for all frames.
%
if size(pupil_mask,3) ~= 1
  error('The pupil mask needs to be the same for all frames.')
end

switch method
  case 'direct'
    P = spdiags(pupil_mask(:), 0, n_ap*n_ap, n_ap*n_ap);
    A = [P*Dx; P*Dy; alpha*speye(n_ap*n_ap)];
    b = [reshape(phase_gradx,n_ap*n_ap,n_frames); reshape(phase_grady,n_ap*n_ap, n_frames); zeros(n_ap*n_ap,n_frames)];
    phase_vecs = A\b;
    pupil_mask = pupil_mask(:,:,ones(1,n_frames));
    phase = pupil_mask .* reshape(phase_vecs, n_ap, n_ap, n_frames);
  otherwise
    error('this only uses a sparse direct factorization for the solve')
end
phase = phase(pad_size+1:m_ap-pad_size, pad_size+1:n_ap-pad_size, :);
