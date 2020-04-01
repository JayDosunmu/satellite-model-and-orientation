function R = SubsampleMatrix(n, ssp)
%
%   R = SubsampleMatrix(n, ssp);
%   
% This function constructs a sparse matrix that implements subsampling
% on an image vector. 
%
% Input:
%   n   - dimension of high resolution grid (that is
%         the high resolution grid is n-by-n pixels)
%         If n is a scalar, it is assumed the high res image is n-by-n
%         If n is a 1-by-2 vector, then it is assumed that the high res
%         image is n(1)-by-n(2)
%   ssp - subsampling parameter
%         If ssp is a scalar, then it is assumed that the subsampling
%         is the same in both x and y directions.
%         If ssp is a 1-by-2 vector, then it is assumed that
%          ssp(1) = subsampling in y-direction (down rows of image)
%          ssp(2) = subsampling in x-direction (across cols of image)
%
% Output: 
%   R   - sparse matrix 
%

%
%  J. Nagy
%  August, 2012
%
%  Modified November, 2015 by J. Nagy
%  This now allows non-integer ssp, and the ssp can be different in
%  the x and y directions.
%

%
%  
if length(n) == 1
    nrows = n; ncols = n;
else
    nrows = n(1); ncols = n(2);
end
if length(ssp) == 1
    ssp_row = ssp; ssp_col = ssp;
else
    ssp_row = ssp(1); ssp_col = ssp(2);
end

nrows_low = round(nrows/ssp_row); % should be an integer, but round just to be sure
ncols_low = round(ncols/ssp_col); % should be an integer, but round just to be sure

%
%  We create a matrices that will do the subsampling of the high res image
%  to obtain the low res images.  
%  First we create the subsampling in the y-direction (subsampling of rows)
%
Drow = sparse(nrows_low, nrows);
idx = round(linspace(1, nrows, nrows_low));
Drow(:,idx) = speye(nrows_low);
%
% The below commented out code is another way to do the subsampling.
%
%mid = round(nrows/2);
%idx1 = round(linspace(1,mid,fix(nrows_low/2)));
%h = min(diff(idx1));
%idx2 = round(linspace(idx1(end)+h,nrows,round(nrows_low/2)));
%Drow(:,[idx1,idx2]) = speye(nrows_low);

%
%  Now create the subsampling in the x-direction (subsampling of cols)
%
Dcol = sparse(ncols_low, ncols);
idx = round(linspace(1, ncols, ncols_low));
Dcol(:,idx) = speye(ncols_low);
%
% The below commented out code is another way to do the subsampling.
%
%mid = round(ncols/2);
%idx1 = round(linspace(1,mid,fix(ncols_low/2)));
%h = min(diff(idx1));
%idx2 = round(linspace(idx1(end)+h,ncols,round(ncols_low/2)));
%Dcol(:,[idx1,idx2]) = speye(ncols_low);

%
%  Before subsampling the high res image, we first apply a low-pass
%  (smoothing) filter over the image.  The next two lines of code
%  create windowed sinc filters.
%
% original is below
% [hrow, ~] = SincKernel(nrows, nrows, 0.1, max(1/ssp_row,0.25));
% [hcol, ~] = SincKernel(ncols, ncols, 0.1, max(1/ssp_col,0.25));

% Changed by Doug to limit the amount of smoothing
% on Apr 8 when working with SOR Ceiling fan data - used 5 for independent
% phases
[hrow, ~] = SincKernel(nrows, nrows, 0.10, max(1/ssp_row,0.25));
[hcol, ~] = SincKernel(ncols, ncols, 0.10, max(1/ssp_col,0.25));



%
%  The SincKernel requires a couple of parameters that change the
%  shape of the kernel.  I'm not sure precisely how these should be 
%  chosen.  Below are some experiments I did, but it seems that
%  the amount of filtering should depend on the subsampling.  If
%  ssp = 1, the high res and low res images should be the same.
%  
%[hrow, Hrow] = SincKernel(nrows, nrows, 0.1, 0.08);
%[hcol, Hcol] = SincKernel(ncols, ncols, 0.1, 0.08);
%
%[hrow, Hrow] = SincKernel(nrows, nrows, 0.1, 0.3);
%[hcol, Hcol] = SincKernel(ncols, ncols, 0.1, 0.3);
%
%[hrow, Hrow] = SincKernel(nrows, nrows, 1, 2);
%[hcol, Hcol] = SincKernel(ncols, ncols, 1, 2);
%
%  The second output of SincKernel implements the filtering with zero
%  boundary conditions.  This seemed to cause slight artifacts on the 
%  boundaries, so we now implement refelctive boundary conditions.
%  The next lines of code do this:
%
crow = (length(hrow)+1)/2;
ccol = (length(hcol)+1)/2;
hrow = [hrow; zeros(nrows-length(hrow),1)];
hcol = [hcol; zeros(nrows-length(hcol),1)];
Arow = sparse(buildToep(hrow, crow) + buildHank(hrow, crow));
Acol = sparse(buildToep(hcol, ccol) + buildHank(hcol, ccol));

%
%  Now that we have the filters and subsampling, we can create one big
%  sparse matrix that performs the transformation from high res image to
%  low res. (The commented out line uses zero boundary conditions.)
%
%R = kron(Dcol*Hcol,Drow*Hrow);
R = kron(Dcol*Acol,Drow*Arow);

%
%--------------------------------------------------
%  This version needs a few subfunctions.
%
function [h, H] = SincKernel(nrows, ncols, fc, b)
%
%  This constructs a windowed sinc filter.  This implementation is 
%  modled after a python code at http://tomroelandts.com
%  
N = ceil(4 / b);
if N/2 == fix(N/2)
    N = N + 1;      
end
n = (0:N-1)';
% Compute sinc filter.
h = p_Sinc(2 * fc * (n - (N - 1) / 2));
% Compute Blackman window.
w = 0.42 - 0.5 * cos(2 * pi * n / (N - 1)) + 0.08 * cos(4 * pi * n / (N - 1));
% Multiply sinc filter with window.
h = h .* w;
% Normalize to get unity gain.
h = h / sum(h);
B = repmat(h.', [min(nrows,ncols),1]);
d = -(N-1)/2:(N-1)/2;
H = spdiags(B, d, nrows, ncols);

%
%-------------------------------------------
%
function T = buildToep(c, k)
%
%  Build a banded Toeplitz matrix from a central column and an index
%  denoting the central column.
%
n = length(c);
col = zeros(n,1);
row = col';
col(1:n-k+1,1) = c(k:n);
row(1,1:k) = c(k:-1:1)';
T = toeplitz(col, row);

%
%-------------------------------------------
%
function H = buildHank(c, k)
%
%  Build a Hankel matrix for separable PSF and reflective boundary
%  conditions.
%
n = length(c);
col = zeros(n,1);
col(1:n-k) = c(k+1:n);
row = zeros(n,1);
row(n-k+2:n) = c(1:k-1);
H = hankel(col, row);