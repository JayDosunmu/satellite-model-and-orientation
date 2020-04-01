function R = SubsampleMatrix(n, ssp)
%
%   R = SubsampleMatrix(n, ssp);
%   
% This function constructs a sparse matrix that implements subsampling
% on an image vector. 
%
% Input:
%   n   - dimension of high resolution grid, assumed square (that is
%         the high resolution grid is n-by-n pixels)
%   ssp - subsampling parameter, which must divide n with no remainder
%
% Output: 
%   R   - sparse matrix 
%

%
%  J. Nagy
%  August, 2012
%

if rem(n,ssp)~=0
  error('n must divide subsampling');
end 

D = kron(speye(n/ssp(1)),ones(1,ssp(1)));
R = kron(D,D)/ssp(1)^2;

