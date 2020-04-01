function [object, IterInfo] = p_ReconstructObject_gdnn(PSF, G, object, MaxIter)
%
%  Projected Gradient Descent Iterative Method for multi-frame 
%  deconvolution, subject to a nonnegativity and volume preserving
%  constraints.
%
%  This can be used to solve least squares problem 
%     min sum ||g_k - conv(PSF_k,object)|| 
%  subject to object >= 0 and sum(object(:)) = v,
%  where g_k and PSF_k are data and PSF for k-th data frame, and v is
%  the data volume.
%
%   Input: 
%     PSF     - 3D array containing PSFs for each frame.
%     G       - 3D array containing 2D-FFT of data frames.
%     object  - initial guess at object -- sum of these pixels will
%                    be conserved.
% 
%   Optional Input:
%     MaxIter - Maximum number of iterations.  Default is n_obj/2, where
%               the object is assumed to have n_obj -by- n_obj pixels.
%
%   Output:
%     object  -  solution
%   IterInfo  -  structure containing some information about the iteration
%                Iter     -  actual number of iterations performed
%                Rnrm     -  norm of the residual at each iteration
%                NE_Rnrm  -  norm of the residual at each iteration
%                StopFlag -  integer that indicates reason for stopping
%                            iteration:
%                               1 ==> Rtol satisfied
%                               2 ==> NE_Rtol satisfied
%                               3 ==> MaxIter reached
%
% J. Nagy, December, 2013
%
%  References: [1] H. Engl, M. Hanke, A. Neubauer. "Regularization of 
%                  Inverse Problems", Kluwer, 2000.
%              [2] P.C. Hansen. "Discrete Inverse Problems: Insight and
%                  Algorithms", SIAM, 2010.
%              [3] C. Vogel. "Computational Methods for Inverse Problems",
%                  SIAM, 2002.
%              [4] Y. Saad. "Iterative Methods for Sparse Linear Systems",
%                  2nd Edition, SIAM, 2003
%   

%
%  Set some tolerances to check for convergence.
%
Rtol = 1e-6;               % absolute residual tolerance
Rtol = norm(G(:)) * Rtol;  % relative residual tolerance

NE_Rtol = 1e-6;            % absolute normal equations tolerance

%
%  We need a little information about number of frames and number of 
%  pixels in the object.
%
nframes = size(G, 3);
n_obj =   size(G, 1);
if nargin == 3
  MaxIter = [];
end
if isempty(MaxIter)
  MaxIter = n_obj/2;
end

%
%  Precompute FFTs of all PSFs
%
H = zeros(n_obj, n_obj, nframes);
for k = 1:nframes
  H(:,:,k) = fft2(fftshift(PSF(:,:,k)));
end

%
% We need to compute the data volume.  Here we use:
%
%data_volume = sum(sum(mean(real(ifft2(G)),3)));

% Commenting this out and getting the data volume from the initial object
% which comes from the abs(WF estimate).
% data_volume = max( sum(sum(real(ifft2(G)),1),2) );
data_volume = sum(object(:));

%
%  If we think of the problem as:
%     min || b - Ax || subject to x>= 0
%  then the major computations are to compute A*vector and A'*vector.
%  But in this particular problem, A* and A'* are convolutions, so 
%  everything is done using FFTs.
%
%  Here, trAb = A'*b 
%
trAb = sum(conj(H).*G, 3);
NE_Rtol = norm(trAb(:)) * NE_Rtol;   % relative normal eq. tolerance

Rnrm = zeros(MaxIter+1, 1);
NE_Rnrm = zeros(MaxIter+1, 1);

%r = b - A*x;
%  
% Again, if we think about this problem as the nonnegative least squares
% problem,
%          min || b - Ax || subject to x>= 0
% then here we use x = object, and r = b - Ax.
%
x = object;
r = G - H.*repmat(fft2(x), [1,1,nframes]);

%
% d = A'*r
%
d = sum(conj(H).*r, 3);

Rnrm(1) = norm(r(:));
NE_Rnrm(1) = norm(d(:));

for k = 1:MaxIter
  if Rnrm(k) <= Rtol
    % stop because residual satisfies ||b-A*x||<= Rtol
    StopFlag = 1;
    break
  end
  if NE_Rnrm(k) <= NE_Rtol
    % stop because normal equations residual satisfies ||A'*b-A'*A*x||<= NE_Rtol
    StopFlag = 2;
    break
  end
  %w = A*d;
  w = H.*repmat(d, [1,1,nframes]);
  % use steepest descent step length
  tau = (d(:)'*d(:))/(w(:)'*w(:));
  %
  %  Do simple line search with backtracking
  %
  di = real(ifft2(d));
  for lstep = 1:10
    xnew1 = x + tau*di;
    xnew = p_gdnn_project(xnew1, data_volume);
    rnew = r - tau*w;
    xdiff = xnew - x;
    if (rnew(:)'*rnew(:)-r(:)'*r(:)) > -(1e-4)*(xdiff(:)'*xdiff(:))/tau;
      tau = tau/2;
    else
      break
    end
  end
  x = xnew;
  r = rnew;
  %d = A'*r;
  d = sum(conj(H).*r, 3);
  Rnrm(k+1) = norm(r(:));
  NE_Rnrm(k+1) = norm(d(:));
end
object = x;

if k == MaxIter
  % Stop because max number of iterations reached
  StopFlag = 3;
else
  k = k - 1;
end
if nargout==2
  IterInfo.Iter = k;
  IterInfo.Rnrm = Rnrm(1:k+1);
  IterInfo.NE_Rnrm = NE_Rnrm(1:k+1);
  IterInfo.StopFlag = StopFlag;
end

