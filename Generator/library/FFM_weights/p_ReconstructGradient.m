function [phaseC_grad, iters, relres, resvec ] = p_ReconstructGradient(A, W, R, phase_grad_meas, n_layers, n_comp, alpha, rtol, Weights)
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

% per Jim's advice to limit the number of iterationss
% MaxIter = 10*(n_frames-1)+1;

alpha = alpha_scale*alpha;

% Weights to mask out bad sub-apertures
if ~isempty(Weights)
    
    % dynamic sub-aperture masking (new version Aug 5, 2018)
    if size(Weights,3)>1
        v_ww = Weights(:,:,1);   WeightMatrix = spdiags(v_ww(:),0,numel(v_ww(:)),numel(v_ww(:)));
        WW = reshape(Weights,[ size(Weights,1)*size(Weights,2) size(Weights,3) ]);
        radv=size(R,1);cadv = size((WeightMatrix*R)*W,2);
        r=1;c=1;
        for j=1:size(Weights,3)
            WW_j = spdiags(WW(:,j),0,numel(WW(:,j)),numel(WW(:,j)));
            Wq(r:r+(radv-1),c:c+(cadv-1)) = (WW_j*R)*W;
            r=r+radv;c=c+cadv;
        end
        A = [ Wq*A  ; alpha*speye(size(A,2)) ]; %= [kron(speye(n_frames),(WeightMatrix*R)*W)*A; alpha*speye(size(A,2))];
        b = [kron(speye(n_frames),WeightMatrix)*phase_grad_meas(:); zeros(size(A,2),1)];
%         r=1;c=1;
%         for j=1:size(Weights,3)
%             WW_j = spdiags(WW(:,j),0,numel(WW(:,j)),numel(WW(:,j)));
%             meas_j = phase_
%             Wq(r:r+(radv-1),c:c+(cadv-1)) = (WW_j*R)*W;
%             r=r+radv;c=c+cadv;
%         end        
    else
        % FFM subap-mask static for entire FFM - Aug 4, 2018
        WeightMatrix = spdiags(Weights(:),0,numel(Weights(:)),numel(Weights(:)));
        A = [kron(speye(n_frames),(WeightMatrix*R)*W)*A; alpha*speye(size(A,2))];
        b = [kron(speye(n_frames),WeightMatrix)*phase_grad_meas(:); zeros(size(A,2),1)];
    end    
% Original FFM     
%      b = [kron(speye(n_frames),WeightMatrix))*phase_grad_meas(:); zeros(size(A,2),1)];
%     b = [WeightMatrix*phase_grad_meas(:); zeros(size(A,2),1)];
else
    A = [kron(speye(n_frames),R*W)*A; alpha*speye(size(A,2))];
    b = [phase_grad_meas(:); zeros(size(A,2),1)];
end
[phaseC_grad, flag, relres, iters, resvec] = lsqr(A, b, rtol, MaxIter);  

%
%  Check this -- looks like we shouldn't need to reshape this into a
%  3D array, because we are using it as a vector.
%
phaseC_grad = reshape(phaseC_grad, n_comp, n_comp, n_layers);

