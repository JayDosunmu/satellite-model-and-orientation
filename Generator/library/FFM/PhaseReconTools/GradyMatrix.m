function Dy = GradyMatrix(n, scale)
%
%  Dy = GradyMatrix(n, scale);
%
% Generates a sparse matrix that models computation of y-gradients
% obtained by a wavefront sensor with Fried geometry.
%
% Important Remarks:
%    1. The geometry is strange -- this computes what really seems like
%       a x-gradient.  However, this code constructs a matrix that is
%       consistent with what is found in the adaptive optics literature.
%       See ...
%    2. However, we do make one change with previously published papers:
%       Our Dy is -1 times what is found in the above cited literature.  
%    3. We are currently implementing this with a reflective boundary
%       condition.
%
% Input:
%   n  -    number of pixels across the image (assumed square)
%
% Optional Input:
%   scale - allows for a normalization scale for specific geometries.
%           Default is to use scale = 1;
%
% Output:
%   Dy -    n^2 x n^2 sparse matrix
%

if nargin == 1
  scale = [];
end
if isempty(scale)
  scale = 1;
end

e = ones(n, 1);
d1 = [e(1:n-1); 2];
d2 = [e(1:n-1); 0];

H = spdiags([-d2, e], 0:1, n, n);
F = .5*spdiags([d1, e], 0:1, n, n);

Dy = scale*kron(H, F);
