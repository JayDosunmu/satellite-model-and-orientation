function phase_layers = ReconstructPhase_Layers(A, W, FFM_input, FFM_output, FFM_options)
%
%    phase_layers = ReconstructPhaseC(FFM_input, FFM_output, FFM_options);
%
%  This function reconstructs the COMPOSITE phases on each layer from 
%  the composite x and y-gradients on each layer. 
%  We assume a Fried geometry of the wavefront sensor.
%
%  Input:
%    FFM_input    - structure that defines input parameters and data 
%                   for the FFM codes
%    FFM_output   - structure that results from running the FFM codes.
%                   In particular, this contains the reconstructed 
%                   composite (sausage region) gradients on each layer.
%
%  Optional Input:
%      FFM_options - structure containing some optional input paramteters.
%                    See FFM.m for more information.
%
%  Output:
%    phase_layers - reconstructed phases on each layer.  This is a 4-D
%                   array, where
%                       phase_layers(:,:,k,L)
%                   contains the phase on layer L for frame k.
%

%
%  J. Nagy
%  December, 2012
%  Modified, July 2015
%  The previous version did not guarantee the frozen flow held on
%  the reconstructed phases for each layer.  This version should
%  fix this problem.
%
%  Modified November, 2015, to allow non square subapertures. (J. Nagy) 
%  So now, ssp can be a scalar (equal diameters in both directions)
%  or a 1-by-2 array, [ssp_rows, ssp_cols]
%
%  Modified April, 2016 to reduce memory requirements.  There are now
%  two new inputs (A and W) which are needed in other parts of the code.
%  So it makes no sense to rebuild them here.

%
%  These first statements just do a quick check to make sure enough
%  input terms have been given; FFM_options is optional.
%  Check options, and set default parameters.
%
%narginchk(2,3)
if nargin == 4
  % Use all defaults for options
  phase_scale = [];
  reg_par_PR = [];
  method_PR = [];
else
  if isfield(FFM_options, 'phase_scale')
    phase_scale = FFM_options.phase_scale;
  else
    phase_scale = [];
  end
  if isfield(FFM_options, 'reg_par_PR')
    reg_par_PR = FFM_options.reg_par_PR;
  else
    reg_par_PR = [];
  end
  if isfield(FFM_options, 'method_PR')
    method_PR = FFM_options.method_PR;
  else
    method_PR = [];
  end
end
if isempty(phase_scale), phase_scale = 1; end
if isempty(reg_par_PR), reg_par_PR = 1e-6; end
if isempty(method_PR), method_PR = 'direct'; end

n_ap = FFM_input.n_ap;
n_frames = FFM_input.n_frames;
n_layers = FFM_input.n_layers;
phaseC_gradx = FFM_output.phaseC_gradx;
phaseC_grady = FFM_output.phaseC_grady;
pupil_mask = FFM_input.pupil_mask;
wind_vecs = FFM_input.wind_vecs;
ssp = FFM_input.ssp;
%
%  Modified November, 2015, to allow non square subapertures.
%  J. Nagy
%
if length(ssp) == 1
    ssp = [ssp, ssp];
end
%
%  Convert wind_vecs into pixel shift values:
%
deltax = ssp(1)*wind_vecs(:,1) .* cos(wind_vecs(:,2));
deltay = ssp(2)*wind_vecs(:,1) .* sin(wind_vecs(:,2));

%
%  We need to find a composite grid size, and the amount of padding needed
%  to get from n_ap to n_conmp.  
%
[n_comp, n_comp_pad] = GetCompositeGridSize(n_ap, n_frames, deltax, deltay);

%
%  In this version, we reconstruct the phase on the composite grid:
%
%pupil_maskC = GetCompositeMask(pupil_mask, wind_vecs, ssp, n_ap, n_comp, n_comp_pad, n_frames, n_layers);
%
%  April 2016 modification: Changed GetCompositeMask.m so that it does not
%  build another big matrix.  Instead it uses the given motion matrix.
%
%
%pupil_maskC = GetCompositeMask(pupil_mask, n_comp, n_comp_pad, n_frames, n_layers);

pupil_maskC = GetCompositeMask(pupil_mask, wind_vecs, ssp(1), n_ap, n_comp, n_comp_pad, n_frames, n_layers);

phaseC_layers = zeros(size(phaseC_gradx));
for k = 1:n_layers
    phaseC_layers(:,:,k) = p_ReconstructPhase(phaseC_gradx(:,:,k), phaseC_grady(:,:,k), pupil_maskC(:,:,k), phase_scale, reg_par_PR);
end


%
%  We need the following matrices to punch out the phases from the
%  composite grid.
%
%  April 2016 modification: These matrices are created in the FFM code,
%  so instead of rebuilding them, pass them in as input.
%
%W = WindowMatrix(n_ap, n_comp, pupil_mask, n_comp_pad+1, n_comp_pad+1);
%A = MotionMatrix(deltax, deltay, n_comp, n_ap, n_frames);

%
%  To punch out, I zero out all other layers, so that the punch out
%  steps adds zeros.  There might be a better way to do this, but it
%  will require a significant rewrite of the code.
%
phase_layers = zeros(n_ap, n_ap, n_frames, n_layers);
for k = 1:n_layers
  phaseC_layers1 = zeros(size(phaseC_layers));
  phaseC_layers1(:,:,k) = phaseC_layers(:,:,k);
  phase_layers(:,:,:,k) = PunchOutFrames(phaseC_layers1, A, W, n_frames, n_layers, n_comp, n_ap);
end
