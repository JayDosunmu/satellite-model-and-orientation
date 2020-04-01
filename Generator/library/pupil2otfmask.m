function mask = pupil2otfmask( pa )

a=ifft2(pa);h=abs(a).^2;H=abs(fft2(h)).^2;
idx = H>1e-5;
mask = zeros(size(pa));
mask(idx)=1.0;