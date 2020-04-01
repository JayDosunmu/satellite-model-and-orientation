function Iout = SimpleImresize(Iin, OutSize)
%
%  Resize given image In to size given by OutSize.
%
%  In = input image to be upsampled
%  OutSize = size you want for the output image.  This can be a scalar
%            or a 1-by-2 array, [nrows, ncols].  If the input is a scalar,
%            say n, then [nrows, ncols] = [n, n].
%
%  Iout = output image, of size 
%
%  The input image can be a 3-d array of image frames, in which case the 
%  output has the same number of frames, each frame resized to the 
%  requested dimension.
% 
[nyin, nxin, nzin] = size(Iin);
switch length(OutSize)
    case 1
        nxout = OutSize; nyout = OutSize; nzout = nzin;
    case 2
        nxout = OutSize(2); nyout = OutSize(1); nzout = nzin;
    case 3
        if nzin ~= OutSize(3)
            error('can only resize 2-d slices of a 3-d array')
        else
            nxout = OutSize(2); nyout = OutSize(1); nzout = nzin;
        end
    otherwise
        error('illegal input for OutSize')
end
Iout = zeros(nyout, nxout, nzout);
[X, Y] = meshgrid(1:nxin,1:nyin);
[Xup, Yup] = meshgrid(linspace(1,nxin,nxout), linspace(1,nyin,nyout));
for k = 1:nzout
    Iout(:,:,k) = interp2(X,Y,Iin(:,:,k),Xup,Yup);
end
