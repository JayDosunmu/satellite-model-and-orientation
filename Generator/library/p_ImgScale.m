function imOut = p_ImgScale(imIn, scale)
% 
%  This function scales an image about its central pixel using a sinc 
%  interpolator.  
%
%  Input:
%     imIn - input image
%     scale - scale parameter
%
%  Note that the original version of this had a third optional input
%  that gave an option for "size" of the output array.  But we were
%  always using "same", so the other options were removed.
%
%  Also, let's assume the size of imIn is square, to simplify stuff.
%    

n = size(imIn,1);
n2 = n/2;
n2f = floor(n2);

kvec = (1:n)';
vone = ones(n,1);
scalei = 1/scale;
scalei2 = scalei*scalei;

s = kvec - n2f;
c = scalei*(kvec - n2);
sincy = p_Sinc(pi*(c*vone' - vone*s'));
sincx = p_Sinc(pi*(vone*c' - s*vone'));

imOut = scalei2*sincy*imIn*sincx;

return

imOut = zeros(size(imIn));

s1 = [1:size(imIn,1)]-floor(size(imIn,1)/2);
s2 = [1:size(imIn,2)]-floor(size(imIn,2)/2);
% Make the sinc interpolators
sincy = zeros(size(imOut,1),size(imIn,1));
sincx = zeros(size(imIn,2),size(imOut,2));
for k = 1:size(imOut,1)
   sincy(k,:) = p_Sinc(pi*((k-size(imOut,1)/2)/scale-s1));
end
for k = 1:size(imOut,2)
   sincx(:,k) = p_Sinc(pi*((k-size(imOut,2)/2)/scale-s2));
end

% Do the interpolation, and normalize to keep data volume the same
imOut = sincy*imIn*sincx/scale^2;

