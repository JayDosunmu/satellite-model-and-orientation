function fval = func_object_move(x,fx,ft,mask)

% fval = func_object_move(x,fx,ft)
% 
% x  = parameters
% fx = object 
% ft = truth object 
%

ft=ft./sum(ft(:));
fx=fx./sum(fx(:));
FX=fft2(fx);
mdim = size(FX,1);
icen = mdim/2+1;


fval    = 0;
im      = zeros(mdim,mdim,2);
xx      = reshape(x,[2 1]);

dr=xx(1);
dc=xx(2);
 
for r=1:mdim
    for c=1:mdim
        im(r,c,1) = -(2*pi*(dc*(c-icen))/mdim);
        im(r,c,2) = -(2*pi*(dr*(r-icen))/mdim);  
    end
end
cexp    = ifftshift(complex(cos(im(:,:,1)+im(:,:,2)),sin(im(:,:,1)+im(:,:,2))));

FXp     = FX.*cexp;
fxp      = real(ifft2(FXp));

res     = (fxp-ft).*mask;
fval    = fval + sum(res(:).^2);
