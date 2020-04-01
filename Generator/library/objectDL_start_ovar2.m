function fx = objectDL_start_ovar2( fx, pa  ) 

tflux = sum(fx(:));
if size(pa,3)>1

    % Potentially this is a DAD run so find pupil with greatest extent
    pixel_ct = zeros(size(pa,3),1);
    for k=1:size(pa,3)
        frm1 = pa(:,:,k);
        frm1 = frm1./max(frm1(:));
        pixel_ct(k) = sum(sum(frm1(:)));
    end
    kk = find( pixel_ct == max(pixel_ct(:)));
    pa1 = pa(:,:,kk(1));
else
    pa1 = pa;
end
a = ifft2(pa1);h=abs(a).^2;h=h./sum(h(:));
Hmodsq  = abs(fftshift(fft2(h))).^2;
idx     = find(Hmodsq > 1e-5);
mmask   = zeros(size(Hmodsq));
mmask(idx) = 1.0;
mmask = ifftshift(mmask);


for j=1:100
    fx_nn   = abs(fx);
    FX_NN   = fft2(fx_nn);
    FX_NN   = FX_NN.*mmask;
    fx      = real(ifft2(FX_NN));
end
fx = abs(fx);
fx = fx.*tflux;