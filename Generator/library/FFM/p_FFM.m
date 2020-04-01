function [FFM_output, A, W, R] = p_FFM(FFM_input, FFM_options)
%
%  FFM_output = p_FFM(FFM_input, FFM_options);
%
%  This function is used to run the FFM stuff for DORA ...
%
%  The p is used to denote that this MATLAB code has been stripped of 
%  any unnecessary computations to try to get it ready for the parallel
%  implementation.
%
%  Here we reconstruct the composite gradients on each layer, but then
%  we punch them out to get gradients on each frame.  Then we recosntruct
%  the phase on each frame.
%
%  Input:
%      FFM_input   - structure containing the following parameters:
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
%  Optional Input:
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
%                                         Default: 1e-6
%         rtol_FFM     Double             Relative residual stopping 
%                                         tolerance for phase reconst.
%                                         This is needed if method_FFM is
%                                         lsqr.
%                                         Default: 1e-6

%
%  Output:
%      FFM_output - structure containing the following parameters:
%
%         Name         Type               Purpose
%         -----------  ---------------    --------------------------------
%         phaseC_gradx 3D double array    contains reconstructed composite 
%                                         x-gradients of wavefront phase 
%                                         on each frame (high resolution)
%         phaseC_grady 3D double array    contains reconstructed composite 
%                                         y-gradients of wavefront phase 
%                                         on each frame (high resolution)
%         A, W, R      Sparse matrices    See below -- these do the motion,
%                                         windowing and subsampling
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
%       deltax = number of pixels the atmosphere shifts to right/left
%       deltay = number of pixels the atmosphere shifts up/down
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

%
%  J. Nagy
%  October, 2013

%
%  These first statements just do a quick check to make sure enough
%  input terms have been given; FFM_options is optional.
%  Check options, and set default parameters.
%
% narginchk(1,2)
if nargin == 1
  % Use all defaults for options
  phase_scale = [];
  reg_par_FFM = [];
  reg_par_PR = [];
  rtol_FFM = [];
else
  if isfield(FFM_options, 'phase_scale')
    phase_scale = FFM_options.phase_scale;
  else
    phase_scale = [];
  end
  if isfield(FFM_options, 'reg_par_FFM')
    reg_par_FFM = FFM_options.reg_par_FFM;
  else
    reg_par_FFM = [];
  end
  if isfield(FFM_options, 'reg_par_PR')
    reg_par_PR = FFM_options.reg_par_PR;
  else
    reg_par_PR = [];
  end
  if isfield(FFM_options, 'rtol_FFM')
    rtol_FFM = FFM_options.rtol_FFM;
  else
    rtol_FFM = [];
  end
end
if isempty(phase_scale), phase_scale = 1; end
if isempty(reg_par_FFM), reg_par_FFM = 1e-3; end
if isempty(reg_par_PR), reg_par_PR = 1e-6; end
if isempty(rtol_FFM), rtol_FFM = 1e-6; end

phase_gradx_meas = FFM_input.phase_gradx;
phase_grady_meas = FFM_input.phase_grady;
pupil_mask = FFM_input.pupil_mask;
ssp = FFM_input.ssp;
n_ap = FFM_input.n_ap;
n_subap = FFM_input.n_subap;
n_frames = FFM_input.n_frames;
n_layers = FFM_input.n_layers;
wind_vecs = FFM_input.wind_vecs;

%
%  Check to make sure the input variables have no obvious errors:
%
if n_ap ~= ssp*n_subap
  error('input values should satisfy n_ap = ssp*n_subap')
end
if size(phase_gradx_meas) ~= size(phase_grady_meas)
  error('input phase gradient arrays should be same size')
end
if n_layers ~= size(wind_vecs,1)
  error('wind_vecs information does not match input number of layers')
end
if n_frames ~= size(phase_gradx_meas,3)
  error('number of frames of input phase gradients does not match n_frames')
end

%
%  Convert wind_vecs into pixel shif\t values:
%
deltax = ssp*wind_vecs(:,1) .* cos(wind_vecs(:,2));
deltay = ssp*wind_vecs(:,1) .* sin(wind_vecs(:,2));
%
%  We need to find a composite grid size, and the amount of padding needed
%  to get from n_ap to n_conmp.  
%
[n_comp, n_comp_pad] = GetCompositeGridSize(n_ap, n_frames, deltax, deltay);

%
%  FFM Gradient Reconstruction
%
%     We first need the following matrices:
%
R = SubsampleMatrix(n_ap, ssp);
W = WindowMatrix(n_ap, n_comp, pupil_mask, n_comp_pad+1, n_comp_pad+1);
A = MotionMatrix(deltax, deltay, n_comp, n_ap, n_frames);

%
%  Now reconstruct composite x-gradients of wavefront phase on each 
%  atmospheric layer:
%
[phaseC_gradx, gradx_iters] = p_ReconstructGradient(A, W, R, ...
  phase_gradx_meas, n_layers, n_comp, reg_par_FFM, rtol_FFM);

%
%  and then reconstruct composite y-gradients of wavefront phase on each
%  atmospheric layer:
%
[phaseC_grady, grady_iters] = p_ReconstructGradient(A, W, R, ...
  phase_grady_meas, n_layers, n_comp, reg_par_FFM, rtol_FFM);

FFM_output.phaseC_gradx = phaseC_gradx;
FFM_output.phaseC_grady = phaseC_grady;
FFM_output.deltax = deltax;
FFM_output.deltay = deltay;