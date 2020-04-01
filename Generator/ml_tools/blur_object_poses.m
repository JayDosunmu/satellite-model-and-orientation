function gp = blur_object_poses( g, Dr0 )

grid_sz = size(g,1);
gp = zeros(size(g));

pdim            = grid_sz./2; % Nyquist sample by default
pupil_pix_pri   = pdim;
pupil_pix_sec   = pdim.*0.09375;

% Generate aperture with Nyquist sampling
pupil1 = make_circle_mask(pdim ,pupil_pix_pri/2 );
pupil2 = make_circle_mask(pdim, pupil_pix_sec/2 );
pupil_mask  = pupil1 - pupil2;
pupil_mask  = padarray(pupil_mask,[ (grid_sz-pdim)/2 (grid_sz-pdim)/2 ],0,'both');
pupil_amp   = pupil_mask;


pupil_amp  = pupil_amp./max(pupil_amp(:));
pupil_amp  = (size(g,1).*pupil_amp)./sqrt(sum(sum(pupil_amp.^2)));        

for k=1:size(g,3)
    phase       = kolmogorov(grid_sz, Dr0 );    
    a   = ifft2(pupil_amp.*exp(sqrt(-1).*phase));
    h   = abs(abs(a).^2);
    gp(:,:,k) = real(ifft2(fft2(h).*fft2(g(:,:,k))));
%     figure(9);imagesc(ifftshift(h));pause(0.5);
end

