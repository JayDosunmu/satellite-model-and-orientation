function  [phase_recon, phase_recon2, phaseC_recon] = TestPhaseRecon(FFM_truth, FFM_input, FFM_options)

pupil_mask = FFM_input.pupil_mask;
n_ap = FFM_input.n_ap;
n_subap = FFM_input.n_subap;
ssp = FFM_input.ssp;
n_frames = FFM_input.n_frames;
n_layers = FFM_input.n_layers;
wind_vecs = FFM_input.wind_vecs;

%
%  Convert wind_vecs into pixel shift values:
%
deltax = ssp*wind_vecs(:,1) .* cos(wind_vecs(:,2));
deltay = ssp*wind_vecs(:,1) .* sin(wind_vecs(:,2));

%
%  We need to find a composite grid size, and the amount of padding needed
%  to get from n_ap to n_conmp.  
%
[n_comp, n_comp_pad] = GetCompositeGridSize(n_ap, n_frames, deltax, deltay);

%
%  First reconstruct phase for each frame:
%
phase_recon = PhaseRecon(FFM_truth.phase_gradx, FFM_truth.phase_grady, FFM_input.pupil_mask);

%
%  Now reconstruct phase on each layer:
%
%phaseC_recon = PhaseRecon(FFM_truth.pupil_maskC.*FFM_truth.phaseC_gradx, ...
%  FFM_truth.pupil_maskC.*FFM_truth.phaseC_grady, FFM_truth.pupil_maskC);
phaseC_recon = PhaseRecon(FFM_truth.phaseC_gradx, FFM_truth.phaseC_grady, FFM_truth.pupil_maskC);

%
%  We want to punch these out, and add to get another version of the
%  reconstructed phase for each frame:
%
W = WindowMatrix(n_ap, n_comp, pupil_mask, n_comp_pad+1, n_comp_pad+1);
A = MotionMatrix(deltax, deltay, n_comp, n_ap, n_frames);

N = n_comp*n_comp;
phaseC_recon2_vec = zeros(N*n_layers,1);
for L = 1:n_layers
  phaseC_recon2_vec((L-1)*N+1:L*N,1) = reshape(phaseC_recon(:,:,L),n_comp*n_comp,1);
end
phase_recon2_vec = kron(speye(n_frames),W)*A*phaseC_recon2_vec;
phase_recon2 = reshape(phase_recon2_vec, n_ap, n_ap, n_frames);
