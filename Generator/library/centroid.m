function c = centroid(im)

% c = centroid(im)
%    Computes the centroid of 2D array im with respect to the center
%    of the array

v = zeros(2,size(im,3));
s = size(im(:,:,1));   
[x y] = meshgrid(1:s(2),1:s(1));
for k=1:size(im,3)
   img = im(:,:,k);
   c(1) = sum(sum(x.*img))/sum(img(:))-s(2)/2;
   c(2) = sum(sum(y.*img))/sum(img(:))-s(1)/2;
%    v(1,k) = c(1);
   % v(2,k) = -c(2);
   
end

