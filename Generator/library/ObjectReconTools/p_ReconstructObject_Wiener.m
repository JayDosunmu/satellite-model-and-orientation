function object = p_ReconstructObject_Wiener(PSF, G, reg_par_WF)
%
%  Wiener filter reconstruction method for multi-frame deconvolution.
%  This can be used to solve least squares problem 
%     min sum ||g_k - conv(PSF_k,object)|| + reg_par_WF*||object||
%  where g_k and PSF_k are data and PSF for k-th data frame.
%
%   Input: PSF         - 3D array containing PSFs for each frame.
%          G           - 3D array containing 2D-FFT of data frames.
%          reg_par_WF  - regularization parameter
% 
%   Output:
%          object      -  reconstructed object
%
% J. Nagy, December, 2013
%   

%
%  We need a little information about number of frames and number of 
%  pixels in the object.
%
nframes = size(G, 3);
n_obj =   size(G, 1);

H = zeros(n_obj, n_obj, nframes);
for k = 1:nframes
  H(:,:,k) = fft2(fftshift(PSF(:,:,k)));
end
F = sum(conj(H).*G,3)./( sum(abs(H).^2,3) + reg_par_WF);
object=real(ifft2(F));
