function [phase_gradx, phase_grady] = PhaseTranspSolve_fast(phase, pupil_mask, scale, alpha, method)
%
%    [phase_gradx, phase_grady] = ReconstructPhase_fast(phase, pupil_mask, alpha);
%
%  This function is needed in the MFBD codes.  In particular, we need
%  to apply the transpose of the phase reconstructor to an array, which
%  has dimensions the same as a phase on a single frame.
%  We assume a Fried geometry of the wavefront sensor.
%
%  Input:
%    phase       - array containing x-gradients of phases.  This can be
%                  either a 2D or a 3D array.  In case of 3D, each 2D
%                  slice, phase(:,:,k), can be a 
%                  particular frame of data.
%    pupil_mask  - pupil mask for phases.  If a 2D array, then it is
%                  that each phase has the same pupil mask.
%
%  Optional Input:
%    scale       - allows for a normalization scale for specific geometries.
%                  Default is to use scale = 1;
%    alpha       - Tikhonov regularization parameter. Default is
%                  alpha = sqrt(eps).
%    method      - character string, either 'direct' or 'lsqr'.
%                  direct:    Tikhonov regularized least squares problems
%                             are solved by MATLAB's backslash operator,  
%                             exploiting sparse matrix tools.
%                  iterative: Tikhonov regularized least squares problems
%                             are solved using the iterative method LSQR.
%                  Default is 'direct'.
%
%  Output:
%    phase_gradx - an array in the phase x-gradient space.
%    phase_grady - an array in the phase y-gradient space
%

%
%  J. Nagy
%  August, 2012
%

%
%  Check inputs and set default values.
%
narginchk(2, 5);
switch nargin
  case 2
    scale = []; alpha = []; method = [];
  case 3
    alpha = []; method = [];
  case 4
    method = [];
end

if isempty(scale)
  scale = 1;
end
if isempty(alpha)
  alpha = sqrt(sqrt(eps));
end
if isempty(method)
  method = 'direct';
end
alpha = scale*alpha;

%
%  This can help reduce some boundary artifacts ...
%
pad_size = 1;
phase = padarray(phase, [pad_size, pad_size], 'both');
pupil_mask = padarray(pupil_mask, [pad_size, pad_size], 'both');

[m_ap, n_ap, n_frames] = size(phase);

Dx = GradxMatrix(n_ap, scale);
Dy = GradyMatrix(n_ap, scale);
P = spdiags(reshape(pupil_mask,n_ap*n_ap,1), 0, n_ap*n_ap, n_ap*n_ap);
A = Dx'*P'*P*Dx + Dy'*P'*P*Dy + alpha*alpha*speye(n_ap*n_ap);

w = A\reshape(phase, n_ap*n_ap, n_frames);
phase_gradx = reshape(P*Dx*w, n_ap, n_ap, n_frames);
phase_grady = reshape(P*Dy*w, n_ap, n_ap, n_frames);

phase_gradx = phase_gradx(pad_size+1:m_ap-pad_size, pad_size+1:n_ap-pad_size, :);
phase_grady = phase_grady(pad_size+1:m_ap-pad_size, pad_size+1:n_ap-pad_size, :);


