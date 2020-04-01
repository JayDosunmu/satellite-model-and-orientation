function phase = PunchOutFrames(phaseC, A, W, n_frames, n_layers, n_comp, n_ap, WA)
%
%  phase = PunchOutFrames(phaseC, A, W, n_frames, n_layers, n_comp, n_ap);
%
%  Input:
%    phaseC   - composite information (could be phases or gradients)
%               for each astmospheric layer
%    A, W     - sparse matrices that model the FFM process of WFS data 
%               collection.
%    n_frames - number of frames of data
%    n_layers - number of atmospheric layers
%    n_comp   - size of composite grid
%    n_ap     - size of aperture
%
%  Output
%    phase    - punnched out versions of input on each frame
%

%
%  J. Nagy
%  August, 2012
%
% modified to precompute the matirx product AW - D.Hope 11/13/2017
%
N = n_comp*n_comp;
phaseC_vec = zeros(N*n_layers,1);
for L = 1:n_layers
  phaseC_vec((L-1)*N+1:L*N,1) = reshape(phaseC(:,:,L),n_comp*n_comp,1);
end
if nargin == 7
     phase_vec = kron(speye(n_frames),W)*A*phaseC_vec; 
else
    phase_vec = WA*phaseC_vec;
end
phase = reshape(phase_vec, n_ap, n_ap, n_frames);