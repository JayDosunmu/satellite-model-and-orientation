function results = image_rmse(fx,ft,center_image,mask_fourier,mask_image)

%% results = image_rmse(estimate,truth,center_image,mask_fourier)

% computes the RMSE in the image/Fourier domains when the truth object is
% known
%
% >> results = image_rmse(estimate,truth)
%
% results = 
% 
%           rmse_fourier: [51,512 double]
%             rmse_image: [512x512 double]
%     rmse_fourier_total: 0.3649
%       rmse_image_total: 0.3649
%     
%

if nargin == 4
    mask = ones(size(ft));
end
      
mdim = size(ft,1);
icen = mdim/2+1;
if nargin < 4 
    mask_fourier = ones(mdim);
end

if center_image 

    % line search to find image shift required to align with truth 
    FX =fft2(fx);    
    options = optimset;options = optimset(options,'Display','iter');
    options = optimset;options = optimset(options,'Display','none','TolFun',1e-6,'TolX',1e-6);    
    x0=rand(2,1);
    [xx,fval,exitflag,output]=fminsearch(@func_object_move,x0(:),options,fx,ft,mask_image);
    dr=xx(1);
    dc=xx(2);

    for r=1:mdim
        for c=1:mdim
            im(r,c,1) = -(2*pi*(dc*(c-icen))/mdim);
            im(r,c,2) = -(2*pi*(dr*(r-icen))/mdim);  
        end
    end
    cexp    = ifftshift(complex(cos(im(:,:,1)+im(:,:,2)),sin(im(:,:,1)+im(:,:,2))));
else   
    FX =fft2(fx);
    cexp = ones(mdim);
end
FXp      = FX.*cexp;
fxp      = real(ifft2(FXp));

FT = fft2(ft);
results.rmse_fourier = mask_fourier.*fftshift(sqrt( abs(FXp - FT).^2./abs(FT).^2));

results.dr = dr;
results.dc = dc;
results.mean_rmse_fourier = sum(results.rmse_fourier(find(mask_fourier==1)))./length(find(mask_fourier==1));
results.rmse_image   = (fxp - ft).^2./ft.^2;
results.rmse_fourier_total =  sqrt(sum(sum(abs(FXp - FT).^2))./sum(sum(abs(FT).^2)));
results.rmse_image_total   =  sqrt(sum(sum( (fxp - ft).^2))./sum(sum(ft.^2)));
results.fx_input   = fx;
results.fx_aligned =  fxp;
results.ft = ft;
results.azi_avg = azi_avg(results.rmse_fourier,ones(size(fx)));


