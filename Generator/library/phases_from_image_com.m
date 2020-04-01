function phases = phases_from_image_com( images )

% phases = phases_from_image_com( images )
%    
mdim = size(images,1);
icen = mdim/2 + 1;

for k=1:size(images,3)
    xx = image_com(images(:,:,k));
    dr=xx(1);
    dc=xx(2);

    for r=1:mdim
        for c=1:mdim
            im(r,c,1) = -(2*pi*(dc*(c-icen))/mdim);
            im(r,c,2) = -(2*pi*(dr*(r-icen))/mdim);  
        end
    end
    phases(:,:,k) = im(:,:,1)+im(:,:,2);
end
