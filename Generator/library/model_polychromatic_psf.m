function psfs_poly  = model_polychromatic_psf(phase, pupil_poly, uber_cfg )


if uber_cfg.dlambda>0
    pupil_max = (1 + uber_cfg.dlambda/(2*uber_cfg.lambda))*uber_cfg.pdim;
    pupil_max = 2*round(pupil_max/2);

    pupil_min = (1 - uber_cfg.dlambda/(2*uber_cfg.lambda))*uber_cfg.pdim;
    pupil_min = 2*round(pupil_min/2);

    %pupild = pupil_min:2:pupil_max;uber
    pupild = pupil_min:pupil_max;

    idx0 = find(pupild == uber_cfg.pdim);

    idx = round(linspace(1, length(pupild), uber_cfg.color_planes));
    if ~any(idx(:) == idx0)
      %  In this case, the pupil diameter spacing strange -- trying to fix
      %  I'm not sure if this conditional statement is needede.  Might be
      %  good to remove.
    %  idx(round(length(idx)/2)) = idx0;
    end
    pupild = pupild(idx);
else     
    pupild = uber_cfg.pdim; 
end



h = zeros(uber_cfg.mdim,uber_cfg.mdim, size(phase,3), length(pupild));
for k = 1:length(pupild)
  poly_scale = uber_cfg.pdim/pupild(k);
  for j = 1:size(phase,3)
    phi  = phase(:,:,j)./poly_scale;
    P    = pupil_poly.*exp(sqrt(-1)*phi);
    psf1 = ifftshift(abs(ifft2(P)).^2);
    %
    % The next function does an interpolation.  Since we always use
    % 'same' in the original code, I've changed this to eliminate the
    % input variable 'same'.
    %
    % psf1 = ImgScale(psf1, poly_scale, 'same');
    psf1 = p_ImgScale(psf1, poly_scale);
    h(:,:,j,k) = psf1;
  end        
end
psfs_poly = squeeze(mean(h,4));
for j=1:size(psfs_poly,3)
    frm = psfs_poly(:,:,j);frm = ifftshift(frm);
    psfs_poly(:,:,j) = frm./sum(frm(:));
end
