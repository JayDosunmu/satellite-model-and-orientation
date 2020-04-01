function [FFM_input, FFM_truth, FFM_options, A, W, R] = TestProblem(n_frames, wind_vecs, noise_level)
%
%  This function will generate some data for a test problem for DORA.
%
%  Optional Input:
%      n_frames     - number of frames used for the simulation.
%                     Default is n_frames = 20
%      wind_vecs    - specifies angle (theta) and magnitude (r) of wind 
%                     velocity for each layer,
%                         wind_vecs(L,1) = r     (L = 1, 2, ..., n_layers)
%                         wind_vecs(L,2) = theta (L = 1, 2, ..., n_layers)
%                     Default is
%                         wind_vecs = [0.05 pi; 0.25 pi/2; 0.15 3*pi/4];
%                     Note that the simulation currently works of at
%                     most three layers.
%      noise_level  - scalar, 0 <= noise_level < 1, which is used to 
%                     determine the percentage of Gaussian white noise
%                     added to the simulated measured data.
%                     Default is noise_level = 0.01 (i.e., 1% noise)
%
%  Output:
%      DORA__ffm_input   - structure containing the following parameters:
%
%         Name         Type               Purpose
%         -----------  ---------------    --------------------------------
%         phase_gradx  3D double array    contains measured x-gradients of 
%                                         wavefront phase on each frame 
%                                         (low resolution)
%         phase_grady  3D double array    contains measured y-gradients of 
%                                         wavefront phase on each frame 
%                                         (low resolution)
%         wind_vecs    2D double array    specifies angle (theta) and 
%                                         magnitude (r) of wind 
%                                         velocity for each layer:
%                                             wind_vecs(L,1) = r   
%                                             wind_vecs(L,2) = theta 
%                                         L = 1, 2, ..., n_layers.
%         pupil_mask   2D logical array   pupil mask for aperture
%         n_subap      Integer            number of pixels across the 
%                                         subaperture
%         n_ap         Integer            number of pixels across the 
%                                         aperture
%         ssp          Integer            subsampling parameter, which 
%                                         should satisfy
%                                              n_ap = ssp*n_subap
%         n_frames     Integer            number of frames of data
%         n_layers     Integer            number of atmospheric layers
%       
%      FFM_truth   - structure containing the following parameters:
%
%         Name         Type               Purpose
%         -----------  ---------------    --------------------------------
%         phase        3D double array    contains true wavefront phases
%                                         for each frame
%         phaseC       3D double array    contains true composite 
%                                         wavefront phases on each 
%                                         atmospheric layer
%         phase_gradx  3D double array    contains true x-gradients of 
%                                         wavefront phase on each frame 
%                                         (high resolution)
%         phase_grady  3D double array    contains true y-gradients of 
%                                         wavefront phase on each frame 
%                                         (high resolution)
%         phaseC_gradx 3D double array    contains true composite 
%                                         x-gradients of wavefront phase 
%                                         on each frame (high resolution)
%         phaseC_grady 3D double array    contains true composite 
%                                         y-gradients of wavefront phase 
%                                         on each frame (high resolution)
%         pupil_maskC  3D logical array   contains pupil mask for
%                                         composite information (phases
%                                         and gradients) on leach layer
%          
%      FFM_options - structure containing the following parameters:
%
%         Name         Type               Purpose
%         -----------  ---------------    --------------------------------
%         phase_scale  Double             parameter used to scale the 
%                                         derivative opertation when 
%                                         reconstructing phase from the
%                                         gradients.
%                                         Default: 1 (i.e., no scaling)
%         reg_par_FFM  Double             Tikhonov regularization parameter 
%                                         for frozen flow model high res
%                                         gradient reconstruction.
%                                         Default: 1e-3
%         reg_par_PR   Double             Tikhonov regularization parameter                              
%                                         for phase reconstruction.
%                                         Default: sqrt(eps)
%         method_FFM   String             Method used to solve the Tikhonov
%                                         regularized least squares problem
%                                         for the FFM gradient
%                                         reconstruction.
%                                         Currently we have only one
%                                         method, 'lsqr', but we might want
%                                         to put in other methods later.
%                                         Default: 'lsqr'
%         method_PR    String             Method used to solve the Tikhonov
%                                         regularized least squares problem
%                                         for the phase reconstruction
%                                         reconstruction. Here we can use
%                                         either 
%                                          'direct': MATLAB's backslash
%                                                    operator, which uses
%                                                    efficient sparse 
%                                                    factorization methods.
%                                          'lsqr'  : iterative LSQR method.
%                                         Default: 'direct'
%         rtol_FFM     Double             Relative residual stopping 
%                                         tolerance for phase reconst.
%                                         This is needed if method_FFM is
%                                         lsqr.
%                                         Default: 1e-6
%         rtol_PR      Double             Relative residual stopping 
%                                         tolerance for phase reconst.
%                                         This is needed if method_PR is
%                                         lsqr.
%                                         Default: 1e-6
%

% Some remarks on notation:
%   * In the codes, I will try to use:
%       k - integer to denote a specific frame of data, k = 1:n_frames
%       L - integer to denote a specific layer, L = 1:n_layers
%       (i,j) - integers denoting a specific pixel in a frame or layer.
%   * In the FFM gradient reconstruction codes, we use pixel shifts
%     given by detlax and deltay.  We assume the information in wind_vecs
%     pertains to the low resolution grid, and so
%       deltax = r*ssp*cos(theta)
%       deltay = r*ssp*sin(theta)
%     where r = wind_vecs(L,1) and theta = wind_vecs(L,2).
%     Note that
%       deltax = number of pixels the atmosphere shifts to the right
%       deltay = number of pixels the atmosphere shifts to the left
%     but this should be on the high resolution grid. So be careful here.
%   * Generally we assume constant speed, but the FFM code does allow for
%     non-constant speed.  Specifically:
%             * Number of rows in deltax (and deltay) = n_layers
%             * If there is only one colum, then it is 
%               assumed that speed is constant.
%             * For non-constant speed, deltax and deltay should
%               have n_frames-1 columns, giving the number
%               of pixels shifting in the x(y)-direction for
%               each frame (the first frame is assumed
%               fixed, so only need 2, 3, ..., n_frames)


switch nargin
  case 0
    n_frames = []; wind_vecs = []; noise_level = [];
  case 1
    wind_vecs = []; noise_level = [];
  case 2
    noise_level = [];
  case 3
    [];
  otherwise
    error('Too many input arguments.')
end

if isempty(n_frames)
  n_frames = 20;
end
if isempty(wind_vecs)
  wind_vecs = [0.25 pi; 0.5 pi/2; 0.65 3*pi/4];
end
if isempty(noise_level)
  noise_level = 0.01;
end

n_layers = size(wind_vecs,1);
  
load phi_comp_true
phaseC_true = phi_comp_true(:,:,1:n_layers)./2;

%
%  We'll assume 32 pixels across subaperture, and ssp = 4
%  (i.e., 128 pixels across aperture), and we'll pick a number of frames:
%
% ssp = 4;
% n_subap = 32;


ssp = 6;
n_subap = 40;

n_ap = ssp*n_subap;

%
%  Convert wind_vecs into pixel shift values:
%
deltax = ssp*wind_vecs(:,1) .* cos(wind_vecs(:,2));
deltay = ssp*wind_vecs(:,1) .* sin(wind_vecs(:,2));

%
%  We need to find a composite grid size, and the amount of padding needed
%  to get from n_ap to n_conmp.  We can then cut down the loaded true
%  composite phases to the largest size that we need for this simulation.
%
[n_comp, n_comp_pad] = GetCompositeGridSize(n_ap, n_frames, deltax, deltay);
phaseC_true = phaseC_true(1:n_comp, 1:n_comp, :);

%
%  With these, we can create pupil masks.  I'm not sure the correct way
%  to do this, but here's what I think:
%    Generally, the pupil mask is an annulus, where the ratio of the 
%    diameters (or radii) of the circles to be 0.25.
%  Our "MakeMask" function uses radius of 1 to generate a circle to the
%  outer edge of the array.  So to geneate the pupil masks, we use:
%
pupil_mask_subap = MakeMask(n_subap, 1, 0.25);
pupil_mask = MakeMask(n_ap, 1, 0.25);
% pupil_mask = ones(n_ap);
%
% Note, for example, if n_ap = 128 (this is the diameter of the outer
% circle in the pupil mask), then the diameter of the inner circle is
% then 0.25*128 = 32 pixels.
% Similarly, if n_subap = 32, then the diameter of the inner circle  is
% then 0.25*32 = 8 pixels.
%
% For another example, in the case of the AEOS data, we used 
% ssp = 6, n_subap = 32, and n_ap = ssp*n_subap = 192.  So in this case,
% the diameter of the inner circle for the pupil mask is 
% 0.25*192 = 48 pixels. 
%

%
%  We need the following matrices:
%
R = SubsampleMatrix(n_ap, ssp);
W = WindowMatrix(n_ap, n_comp, pupil_mask, n_comp_pad+1, n_comp_pad+1);
A = MotionMatrix(deltax, deltay, n_comp, n_ap, n_frames);

%
%  Get punched out versions of the wavefront on each frame
%
N = n_comp*n_comp;
phase_true_comp_vec = zeros(N*n_layers,1);
for L = 1:n_layers
  phase_true_comp_vec((L-1)*N+1:L*N,1) = reshape(phaseC_true(:,:,L),n_comp*n_comp,1);
end

phase_true_vec = kron(speye(n_frames),W)*A*phase_true_comp_vec;
phase_true = reshape(phase_true_vec, n_ap, n_ap, n_frames);
%
% The next step in the simulation is to construct gradients.
%
Dx = GradxMatrix(n_comp);
Dy = GradyMatrix(n_comp);
 
phaseC_gradx_true = zeros(size(phaseC_true));
phaseC_grady_true = zeros(size(phaseC_true));
for L = 1:n_layers
  phaseC_gradx_true(:,:,L) = reshape(Dx*reshape(phaseC_true(:,:,L),n_comp*n_comp,1),n_comp,n_comp);
  phaseC_grady_true(:,:,L) = reshape(Dy*reshape(phaseC_true(:,:,L),n_comp*n_comp,1),n_comp,n_comp);
end

%
% Get punched out versions of these gradients, and the corresponding low
% resolution versions:
%
N = n_comp*n_comp;
phaseC_gradx_true_vec = zeros(N*n_layers,1);
phaseC_grady_true_vec = zeros(N*n_layers,1);
for L = 1:n_layers
  phaseC_gradx_true_vec((L-1)*N+1:L*N,1) = reshape(phaseC_gradx_true(:,:,L),n_comp*n_comp,1);
  phaseC_grady_true_vec((L-1)*N+1:L*N,1) = reshape(phaseC_grady_true(:,:,L),n_comp*n_comp,1);
end
phase_gradx_true_vec = kron(speye(n_frames),W)*A*phaseC_gradx_true_vec;
phase_grady_true_vec = kron(speye(n_frames),W)*A*phaseC_grady_true_vec;
phase_gradx_meas_vec = kron(speye(n_frames),R)*phase_gradx_true_vec;
phase_grady_meas_vec = kron(speye(n_frames),R)*phase_grady_true_vec;

phase_gradx_true = reshape(phase_gradx_true_vec, n_ap, n_ap, n_frames);
phase_grady_true = reshape(phase_grady_true_vec, n_ap, n_ap, n_frames);
phase_gradx_meas = reshape(phase_gradx_meas_vec, n_subap, n_subap, n_frames);
phase_grady_meas = reshape(phase_grady_meas_vec, n_subap, n_subap, n_frames);

%
%  Now add noise to the measured gradients -- using the same noise level,
%  but a different realization for each frame.
%
for k = 1:n_frames
  phase_gradx_meas(:,:,k) = phase_gradx_meas(:,:,k) + ...
    WhiteNoise(phase_gradx_meas(:,:,k), noise_level, k);
  phase_grady_meas(:,:,k) = phase_grady_meas(:,:,k) + ...
    WhiteNoise(phase_grady_meas(:,:,k), noise_level, n_frames+k);
end
%(A, pupil_mask, n_comp, n_comp_pad, n_frames, n_layers)
pupil_maskC = GetCompositeMask(A,pupil_mask, n_comp, n_comp_pad, n_frames, n_layers);

FFM_input.phase_gradx = phase_gradx_meas;
FFM_input.phase_grady = phase_grady_meas;
FFM_input.pupil_mask = pupil_mask;
FFM_input.wind_vecs = wind_vecs;
FFM_input.ssp = ssp;
FFM_input.n_ap = n_ap;
FFM_input.n_subap = n_subap;
FFM_input.n_frames = n_frames;
FFM_input.n_layers = n_layers;
FFM_input.subap_mask = pupil_mask_subap;

FFM_truth.phase = phase_true;
FFM_truth.phaseC = phaseC_true;
FFM_truth.phase_gradx = phase_gradx_true;
FFM_truth.phase_grady = phase_grady_true;
FFM_truth.phaseC_gradx = phaseC_gradx_true;
FFM_truth.phaseC_grady = phaseC_grady_true;
FFM_truth.pupil_maskC = pupil_maskC;

FFM_options.phase_scale = 1;
FFM_options.reg_par_FFM = 1e-3;
FFM_options.reg_par_PR = sqrt(eps);
FFM_options.method_FFM = 'lsqr';
FFM_options.method_PR = 'direct';
FFM_options.rtol_FFM = 1e-6;
FFM_options.rtol_PR = [];



