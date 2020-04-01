function [ out_arr ] = padarray(arr,v_dim,value,cwhat)

% padarray(z_poly(:,:,q),[xdim xdim],0,'both');
%
% 11/4/14 - fixed to pad a cube instead of just a single frame 
%
out_arr=zeros(size(arr,1)+2*v_dim(1),size(arr,2)+2*v_dim(2),size(arr,3));
out_arr(v_dim(1)+1:v_dim(1)+size(arr,1),v_dim(2)+1:v_dim(2)+size(arr,2),:)=arr;



end

