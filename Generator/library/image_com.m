function [vecs,phase0] = image_com(im)

% c = image_com(im)
%    Computes the centroid of 2D array im with respect to the center
%    of the array

v = zeros(2,size(im,3));
s = size(im(:,:,1));   
[x y] = meshgrid(1:s(2),1:s(1));

mdim    = size(im,1);
icen    = mdim/2+1;

%** calculate the centroids 
pix_num  = -mdim/2:1:(mdim-1)/2; 
vecs = zeros(2,size(im,3));

phase0 = zeros(size(im));
for k=1:size(im,3)
    img1_sum = sum(im(:,:,k),1);
    img2_sum = sum(im(:,:,k),2);
    c(1) = sum(pix_num(:).*img1_sum(:))/sum(img1_sum);
    c(2) = sum(pix_num(:).*img2_sum(:))/sum(img2_sum);
    vecs(:,k) = [ c(1) c(2) ]';
    [x,y]=meshgrid(-mdim/2:1:mdim/2-1,-mdim/2:1:mdim/2-1);
    im1 = (2*pi*((vecs(1,k).*x))/mdim);
    im2 = (2*pi*((vecs(2,k).*y))/mdim);
    phase0(:,:,k) = -(im1+im2);
       
end

