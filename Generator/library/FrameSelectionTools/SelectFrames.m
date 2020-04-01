function [selected_idx, g_sorted_all, g_selected, g_power, g_SNR_mask] = SelectFrames(g, p)
%
%  This function will sort a set of image data frames in order of 
%  highest fill in the Fourier plane. The cumulative fill factor of the 
%  Fourier plane is shown against each frame in the ordered list. The 
%  user can decide on the trade off between fill factor and number of 
%  frames.
%
% Input:   g - stack of image frames, stored in a 3D array
%
% Optional Input:
%          p - scalar, 0 <= p <= 1, specifies percent of fill
%              in Fourier plane. Default value is p = 1.
%
% Output:  selected_idx - index of selected frames
%          g_selected   - these are the selected frames, in the order
%                         given by frame_idx.  That is, 
%                            g_selected = g(frame_idx)
%
% Optional Output:
%          g_sorted_all - ordered stack of all input image frames
%          g_power      - power spectrum of g, scaled by noise level
%          g_SNR_mask   - mask used to determine Fourier plane
%
%

% Set default input values, and get some parameters
if nargin == 1
  p = 1;
end
n_frames = size(g,3);
n_row = size(g,1);

%
%  First compute the power spectrum, and fftshift.  Before doing this,
%  though, scale so that the max pixel value of each frame is 1.
%
g_power = zeros(size(g));
for k = 1:n_frames
  tt = g(:,:,k)/max(max(g(:,:,k)));
  g_power(:,:,k) = fftshift(abs(fft2(tt)).^2);
end

%
%  Now we need to estimate the noise level, and get the SNR of each 
%  pixel in each image frame.
%
%  To get an estimate of the noise, use the corners of the image
%  frames.  To get the corners, we use a circular mask and then find
%  the mean of pixels outside the circular mask.
%
%  Then divide each frame by these mean values -- this gives us an
%  SNR for each pixel in each frame.
%
n_row = size(g_power,1);

Corner_Mask = MakeMask(n_row, 1);
Corner_Mask = Corner_Mask(:,:,ones(1,n_frames));
g_power_Masked = ~Corner_Mask.*g_power;
g_SNR = zeros(size(g));
for k = 1:n_frames
  g_power_Masked_k = g_power_Masked(:,:,k);
  NoiseEst_frame_k = mean(g_power_Masked_k(g_power_Masked_k~=0));
  nn(k,1) = NoiseEst_frame_k;
  g_SNR(:,:,k) = g_power(:,:,k)/NoiseEst_frame_k;
end
%
%  Now find mask for the object.  
%
%  For this part, first sum all data frames to get a "mean" SNR image
%  (dividing by the number of image isn't necessary to get a true mean)
%
g_SNR = Corner_Mask.*g_SNR;
g_SNR_sum = sum(g_SNR, 3);
%
%  Next we need to find a threshold for the mask.  We do this by
%  sorting pixels in the mean SNR image, and determine where these
%  values flatten out.  That is, we try to find a point of maximum
%  curvature.  To find this "corner" point, we use a technique very
%  similar to what is done for L-curve regularization (see the m-file
%  FindCorner.m for more information).
%
v = sort(g_SNR_sum(:),'descend');
[x_c, y_c] = FindCorner(1:length(v), log(v+1));

%threshold = exp(y_c-1);
threshold = exp(y_c)-1;
g_SNR_mask = g_SNR_sum;
g_SNR_mask(g_SNR_sum<=threshold) = 0;
g_SNR_mask(g_SNR_mask>0) = 1;
g_SNR_mask_all = g_SNR_mask(:,:,ones(1,n_frames));

g_SNR = g_SNR_mask_all.*g_SNR;

Sort_Mask = zeros(size(g_SNR));

max_SNR = max(g_SNR,[],3);
for k = 1:n_frames
  Sort_Mask(:,:,k) = max_SNR == g_SNR(:,:,k);
end
Sort_Mask = g_SNR_mask_all.*Sort_Mask;

nonzeros = zeros(n_frames,1);
for k = 1:n_frames
  nonzeros(k,1) = sum(sum(Sort_Mask(:,:,k)));
end

%
% Sort by number of nonzeros in Sort_Mask
%
[s, sidx] = sort(nonzeros,'descend');
g_sorted_all = g(:,:,sidx);
Sort_Mask = Sort_Mask(:,:,sidx);

num_pixels = sum(g_SNR_mask(:));
sum_pixels = 0;
for k = 1:n_frames
  sum_pixels = sum_pixels + sum(sum(Sort_Mask(:,:,k)));
  if sum_pixels >= p*num_pixels;
    break
  end
end
g_selected = g_sorted_all(:,:,1:k);
selected_idx = sidx(1:k);