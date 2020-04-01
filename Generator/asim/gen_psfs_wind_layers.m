function gen_psfs_wind_layers(Dr0,v_x,v_y,uber_cfg)

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


Nseq_length = 200;
arr = zeros(uber_cfg.mdim,uber_cfg.mdim,Nseq_length);
kk_ct = 0;
for k1=1:length(v_x)
    
    vmag = sqrt(v_x(k1).^2 + v_y(k1).^2);
    if vmag>0
        tic
        grid_sz = 4*(2^ceil(log2(uber_cfg.pdim*vmag)));
        xscale = grid_sz/uber_cfg.mdim;
        
        % Generate layer-1 
        phase2  = kolmogorov(grid_sz, xscale*Dr0 );
        p0      = (grid_sz - uber_cfg.pdim)/2;
        v_co    = (p0+1):(p0+uber_cfg.pdim);
        cube    = zeros(uber_cfg.pdim,uber_cfg.pdim,Nseq_length);

        % Propagate the large layer using randomly chosen wind vectors 
        tmp     = phase2;
        for k=1:Nseq_length
            cube(:,:,k) = tmp(v_co,v_co);
            if v_x(k1)>0
                tmp=shiftr(tmp,0,abs(v_x(k1)));
            elseif v_x(k1)<0
                tmp=shiftl(tmp,0,abs(v_x(k1)));
            end
            if v_y(k1)>0
                tmp=shiftu(tmp,0,abs(v_y(k1)));
            elseif v_y(k1)<0
                tmp=shiftd(tmp,0,abs(v_y(k1)));
            end
        end
        cube        = cube.*uber_cfg.pupil_mask;
        
        % Generate a second layer with a differnt random veloicty 
        phase       = padarray(cube,[pad_dim pad_dim],0,'both');
        psf         = zeros(size(phase));



        for k=1:size(phase,3)
            a   = ifft2(pupil_amp.*exp(sqrt(-1).*phase(:,:,k)));
            h   = abs(ifftshift(abs(a).^2));
            psf(:,:,k) = h;
        end

        arr = abs(fftshift(fft2(psf))).^2;
        save(['mc' num2str(k1) '_psf_modulus_squared.mat' ],'arr','psf');
        tt=toc;
        kk_ct = kk_ct+1;
        fprintf([ 'MC: ' num2str(kk_ct) ' \t v = (' num2str(v_x(k1)) ',' num2str(v_y(k1)) ') \t elapsed time: ' num2str(tt,'%5.3f seconds') '\n' ]);
    end
end



