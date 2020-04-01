function vecs = sm_centroid(im,T)

% c = sm_centroid(im,T)
%    Computes the centroid of 2D array im with respect to the center
%    of the array

if nargin==1;T=0;end
v = zeros(2,size(im,3));
s = size(im(:,:,1));   
[x y] = meshgrid(1:s(2),1:s(1));

mdim    = size(im,1);
icen    = mdim/2;
mask    = ones(mdim);

%** calculate the centroids 
pix_num  = -mdim/2:1:(mdim-1)/2; 
vecs = zeros(2,size(im,3));
for k=1:size(im,3)
    if T>0
        frm = im(:,:,k);
        idx=find(im(:,:,k)<T.*max(frm(:)));
        mask=zeros(mdim);mask(idx)=1.0;                
    end
    
    img1_sum = sum(im(:,:,k).*mask,1);
    img2_sum = sum(im(:,:,k).*mask,2);
    c(2) = -sum(pix_num.*img1_sum)/sum(img1_sum);
    c(1) = -sum(pix_num*img2_sum)/sum(img2_sum);
    vecs(:,k) = [ c(1) c(2) ]';
    
end

