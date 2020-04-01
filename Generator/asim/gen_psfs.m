function gen_psfs(Dr0,v_x,v_y,uber_cfg)

% gen_psfs(Dr0,v_x,v_y) 
% 
% Dr0    Turbulence Strength
% v_x 	 x-shift
% v_y    y-shift 
%
pad_dim     = (uber_cfg.mdim - uber_cfg.pdim)/2;

pupil_amp = padarray(uber_cfg.pupil_mask,[ pad_dim pad_dim],0,'both');
pupil_amp = pupil_amp./max(pupil_amp(:));
pupil_amp = (uber_cfg.mdim.*pupil_amp)./sqrt(sum(sum(pupil_amp.^2)));


Nseq_length = 100;
for k1=1:length(v_x)
    
    vmag = sqrt(v_x(k1).^2 + v_y(k1).^2);

    tic
    grid_sz = 2*(2^ceil(log2(uber_cfg.pdim*vmag)));
    phase2  = kolmogorov(grid_sz, Dr0 );
    p0      = (grid_sz - uber_cfg.pdim)/2;
    v_co    = (p0+1):(p0+uber_cfg.pdim);
    cube    = zeros(uber_cfg.pdim,uber_cfg.pdim,Nseq_length);

    tmp     = phase2;
    for k=1:100
        cube(:,:,k) = tmp(v_co,v_co);
        if v_x>0
            tmp=shiftr(tmp,0,abs(v_x));
        elseif v_x<0
            tmp=shiftl(tmp,0,abs(v_x));
        end
        if v_y>0
            tmp=shiftu(tmp,0,abs(v_y));
        elseif v_y<0
            tmp=shiftd(tmp,0,abs(v_y));
        end
    end

    phase       = padarray(cube,[pad_dim pad_dim],0,'both');
    psf         = zeros(size(phase));



    for k=1:size(phase,3)
        a   = ifft2(pupil_amp.*exp(sqrt(-1).*phase(:,:,k)));
        h   = abs(ifftshift(abs(a).^2));
        psf(:,:,k) = h;
    end

    save(['mc' num2str(k1) '_psfs.mat' ],'psf');

    tt=toc;
    fprintf([ 'MC: ' num2str(k1) ' \t v = (' num2str(v_x(k1)) ',' num2str(v_y(k1)) ') \t elapsed time: ' num2str(tt,'%5.3f seconds') '\n' ]);
%     fitswrite(psf,['mc' num2str(k1) '_psfs.fits' ]);
end



