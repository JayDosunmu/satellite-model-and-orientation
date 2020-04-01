function [ out_arr ] = padarray(arr,pad_dim,where,where2)

if  nargin > 3 | strcmpi(where,'both')
    nslices = size(arr,3);
    ndim = pad_dim(1) + size(arr,1) + pad_dim(2);
    ndim2 = size(arr,1);
    out_arr = zeros(ndim,ndim,nslices);
    for kk=1:nslices
        out_arr(pad_dim(1)+1:(pad_dim(1)+1 + ndim2-1),pad_dim(1)+1:(pad_dim(1)+1 + ndim2-1),kk) = arr(:,:,kk);
    end
end
%  phasek_pad = padarray(padarray(phasek, pad_size_pre, 'pre'), pad_size_post, 'post');
nslices = size(arr,3);
if strcmpi(where,'pre')   
        out_arr = zeros( size(arr,1) + pad_dim(1),size(arr,2) + pad_dim(2),nslices);
        for kk=1:nslices
            out_arr(pad_dim+1:end,pad_dim+1:end,kk) = arr(:,:,kk);
        end        
end
if strcmpi(where,'post')
    out_arr = zeros( size(arr,1) + pad_dim(1),size(arr,2) + pad_dim(2),nslices);
    for kk=1:nslices 
        out_arr(1:size(arr,1),1:size(arr,2),kk) = arr(:,:,kk);
    end
end


end

