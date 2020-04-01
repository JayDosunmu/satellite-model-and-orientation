function [phaseC_grad, iters] = p_ReconstructGradient(A, W, R, phase_grad_meas, n_layers, n_comp, alpha, rtol)
%
%  [phaseC_grad, iters] = ReconstructGradient(A, W, R, ...
%         phase_grad_meas, n_layers, n_comp, alpha, method, rtol);
%    
%  Use the frozen flow model (FFM) to reconstruct high resolution 
%  gradients from low resolution measurements. 
%
%  The p is used to denote that this MATLAB code has been stripped of 
%  any unnecessary computations to try to get it ready for the parallel
%  implementation.
%
%  Input:
%    A, W, R          - sparse matrices that model the FFM process of
%                       WFS data collection.
%    phase_grad_meas  - 3D array containing measured (low resolution) 
%                       x or y-gradients of phases, where (:,:,k) is the
%                       k-th frame of data.
%    n_layers         - number of atmospheric layers
%    n_comp           - size of the composite grid
%
%    alpha            - Tikhonov regularization parameter. Default is
%                       alpha = 1e-3.
%    rtol             - residual stopping tolerance for LSQR.
%                       Default is 1e-6.
%
%  Output:
%    phaseC_grad      - 3D array containing reconstructed (high resolution) 
%                       x or y-gradients of phases on composite grid, 
%                       where (:,:,L) is the L-th atmospheric layer
%
%  Optional Output:
%    iters            - number of LSQR iterations needed to find least
%                       squares solution.
%

%
%  J. Nagy
%  October, 2013
%

alpha_scale = sqrt(max(A*ones(size(A,2),1))*max(A'*ones(size(A,1),1)));

n_frames = size(phase_grad_meas, 3);
MaxIter = 100*n_comp;
alpha = alpha_scale*alpha;
A = [kron(speye(n_frames),R*W)*A; alpha*speye(size(A,2))];
b = [phase_grad_meas(:); zeros(size(A,2),1)];
[phaseC_grad, flag, relres, iters, resvec] = lsqr(A, b, rtol, MaxIter);
%
%  Check this -- looks like we shouldn't need to reshape this into a
%  3D array, because we are using it as a vector.
%
phaseC_grad = reshape(phaseC_grad, n_comp, n_comp, n_layers);

