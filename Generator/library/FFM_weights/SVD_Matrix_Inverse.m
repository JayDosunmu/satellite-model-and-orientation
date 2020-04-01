function invM = SVD_Matrix_Inverse(M, conditionNumber)
%
%    invM = SVD_Matrix_Inverse(M, conditionNumber);
%
%  This function inverts the input matrix using SVD with the specified
%  condition number (ratio of highest tolowest singular values retained).
%
%  Input:
%    M - 2D matrix to be inverted
%
%  Optional Input:
%    conditionNumber - ratio of highest tolowest singular values retained
%                      (default: 100.0)
%
%  Output:
%    invM - SVD inverse of input matrix M
%
%  M. Milton - March 30, 2018

%
%  Check inputs and set default values.
%
if nargin < 2
    conditionNumber = 1.0e2;
end

% Compute inverse of M using SVD
[U, S, Vtranspose] = svd(M);
[Srows, Scols] = size(S);
Sdiag = diag(S);
Sdiag((Sdiag(1)./Sdiag)>conditionNumber)=0;
Sdiag=1./Sdiag;
Sdiag(isinf(Sdiag))=0;
Sdiag(isnan(Sdiag))=0;
Sinv = zeros(Scols, Srows);
Sinv(1:Scols,1:Scols)=diag(Sdiag);
invM = Vtranspose * Sinv * transpose(U);
