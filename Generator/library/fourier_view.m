function fourier_view( fname )

ft=fitsread(fname);
cube = zeros(size(ft));
fname2 = [fname(1:findstr(fname,'.fits')-1) '_modsq.fits' ];
for k=1:size(cube,3)
    cube(:,:,k) = abs(fftshift(fft2(ft(:,:,k)))).^2;   
end

fitswrite(cube,fname2)
fprintf(['Fourier modsq -> ' fname2 '\n\n' ]);