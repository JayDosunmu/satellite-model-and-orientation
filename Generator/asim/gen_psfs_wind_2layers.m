function gen_psfs_wind_layers(Dr0,wind_v,uber_cfg)

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


v_x1 = wind_v(:,1,1);v_y1 = wind_v(:,2,1);
v_x2 = wind_v(:,1,2);v_y2 = wind_v(:,2,2);


Nseq_length = 200;
arr = zeros(uber_cfg.mdim,uber_cfg.mdim,Nseq_length);
kk_ct = 0;
for k1=1:size(wind_v,1)
    
    L1_v = [ v_x1(k1) v_y1(k1) ];
    L2_v = [ v_x2(k1) v_y2(k1) ];
    
    vmag1   = sqrt(v_x1(k1).^2 + v_y1(k1).^2);
    vmag2   = sqrt(v_x2(k1).^2 + v_y2(k1).^2);
    vmag    = max(vmag1,vmag2);
    

    if vmag1>0 & vmag2 >0
        tic
        grid_sz = 4*(2^ceil(log2(uber_cfg.pdim*vmag)));
        xscale = grid_sz/uber_cfg.mdim;
        
        % Generate layer-1 
        phase2  = kolmogorov(grid_sz, xscale.*Dr0 );
        p0      = (grid_sz - uber_cfg.pdim)/2;
        v_co    = (p0+1):(p0+uber_cfg.pdim);
        cube1    = zeros(uber_cfg.pdim,uber_cfg.pdim,Nseq_length);cube2 = cube1;

        % Propagate the large phase screen using randomly chosen L1-wind vectors 
        tmp     = phase2;
        for k=1:Nseq_length
            cube1(:,:,k) = tmp(v_co,v_co);
            if v_x1(k1)>0
                tmp=shiftr(tmp,0,abs(v_x1(k1)));
            elseif v_x1(k1)<0
                tmp=shiftl(tmp,0,abs(v_x1(k1)));
            end
            if v_y1(k1)>0
                tmp=shiftu(tmp,0,abs(v_y1(k1)));
            elseif v_y1(k1)<0
                tmp=shiftd(tmp,0,abs(v_y1(k1)));
            end
        end
        cube1        = cube1.*uber_cfg.pupil_mask;
               

        % Generate layer-2 
        grid_sz = 4*(2^ceil(log2(uber_cfg.pdim*vmag)));
        xscale = grid_sz/uber_cfg.mdim;        
        phase2  = kolmogorov(grid_sz, xscale.*Dr0 );
        p0      = (grid_sz - uber_cfg.pdim)/2;
        v_co    = (p0+1):(p0+uber_cfg.pdim);

        % Propagate the large phase screen using randomly chosen L1-wind vectors 
        tmp     = phase2;
        for k=1:Nseq_length
            cube2(:,:,k) = tmp(v_co,v_co);
            if v_x2(k1)>0
                tmp=shiftr(tmp,0,abs(v_x2(k1)));
            elseif v_x2(k1)<0
                tmp=shiftl(tmp,0,abs(v_x2(k1)));
            end
            if v_y2(k1)>0
                tmp=shiftu(tmp,0,abs(v_y2(k1)));
            elseif v_y2(k1)<0
                tmp=shiftd(tmp,0,abs(v_y2(k1)));
            end
        end
        cube2       = cube2.*uber_cfg.pupil_mask;
        
        % Generate a second layer with a differnt random veloicty 
        phase       = padarray(cube1+cube2,[pad_dim pad_dim],0,'both');
        psf         = zeros(size(phase));

        for k=1:size(phase,3)
            a   = ifft2(pupil_amp.*exp(sqrt(-1).*phase(:,:,k)));
            h   = abs(ifftshift(abs(a).^2));
            psf(:,:,k) = h;
        end

        arr = abs(fftshift(fft2(psf))).^2;
        save(['mc' num2str(k1) '_2L_psf_modulus_squared.mat' ],'arr','psf','phase','cube1','cube2','L1_v','L2_v');
        tt=toc;
        kk_ct = kk_ct+1;
        sL1 = [ 'L1_v = (' num2str(v_x1(k1)) ',' num2str(v_y1(k1)) ')' ];
        sL2 = [ 'L2_v = (' num2str(v_x2(k1)) ',' num2str(v_y2(k1)) ')' ];
        
        fprintf([ 'MC: ' num2str(kk_ct) ' \t' sL1 '  ' sL2 '\t elapsed time: ' num2str(tt,'%5.3f seconds') '\n' ]);
    end
end



