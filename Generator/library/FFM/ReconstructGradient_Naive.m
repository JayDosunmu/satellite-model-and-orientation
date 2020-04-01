function [phase_gradx, phase_grady] = ReconstructGradient_Naive(phase_gradx_meas, phase_grady_meas, n_ap, ssp)
%
%  [phase] = ReconstructGradient_Naive(phase_grad_measa);
%    
%  Use up-sampling interpolation to naively reconstruct phase gradients.
%
%  Input:
%    phase_grad_meas  - 3D array containing measured (low resolution) 
%                       x or y-gradients of phases, where (:,:,k) is the
%                       k-th frame of data.
%  Output:
%    phase_gradx      - 3D array containing reconstructed (high resolution) 
%                       x-gradients of phases where (:,:,k) is the 
%                       k-th frame of data.
%    phase_grady      - 3D array containing reconstructed (high resolution) 
%                       y-gradients of phases where (:,:,k) is the 
%                       k-th frame of data.
%

%
%  J. Nagy
%  January, 2014
%

nframes = size(phase_gradx_meas, 3);

phase_gradx = zeros(n_ap, n_ap, nframes);
phase_grady = zeros(n_ap, n_ap, nframes);
for k = 1:nframes
  phase_gradx(:,:,k) = imresize(phase_gradx_meas(:,:,k),ssp);
  phase_grady(:,:,k) = imresize(phase_grady_meas(:,:,k),ssp);
end

