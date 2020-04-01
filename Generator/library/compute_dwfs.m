function [afx,v_fit] = compute_dwfs(SM_params,G,aphase,wfs_transpose,wfs_rot_angle)


mdim        = size(G,1);
Gmodsq      = fftshift(abs(G).^2);
noise       = mean(Gmodsq(1:10,1:10,:));err = mean(noise(:));
SNR         = Gmodsq./err;
idx         = find(SNR>10);
mask        = zeros(size(SNR));
mask(idx)   = 1.0;
for k=1:size(mask,3);mask(:,:,k) = ifftshift(mask(:,:,k));end


%% Vary pupil size and phase scaling to get the best object
v_size =  108; %:2:192; %[60:2:192 ];
v_scale = 1:0.1:4;

v_max = zeros(length(v_scale),length(v_size));v_fit = v_max;
afx   = zeros(mdim,mdim,length(v_scale),length(v_size));
H=zeros(size(G));

v_int = 16;

for ss=1:length(v_size)
    
    xx          = (SM_params.Diam_sec./SM_params.Diam_pri);
    pupil_outer = dora_make_mask(v_size(ss),v_size(ss));
    r0          = round( xx*v_size(ss)/2);
    pupil_inner = dora_make_mask(v_size(ss),r0);
    pupil_mask  = pupil_outer - pupil_inner;

    % Make Pupil amplitude mask 
    pad_size = (mdim - v_size(ss))/2;
    pad_size_pre  = fix([pad_size, pad_size]);
    pad_size_post = pad_size_pre;

    pupil_amp  =  padarray(padarray(pupil_mask, pad_size_pre, 'pre'), pad_size_post, 'post');
    pupil_amp  = pupil_amp./max(pupil_amp(:));
    pupil_amp  = (mdim.*pupil_amp)./sqrt(sum(sum(pupil_amp.^2)));

    for aa=1:length(v_scale)
        for k=1:size(G,3)
            hh = zeros(mdim);
            for jp=1:length(v_int)
                jj          = v_int(jp);
               phi         = aphase(:,:,jj,k);                            
%                 phiR         = rot90(phi,round(wfs_rot_angle/90));  
                       phiR         = imrotate(phi,wfs_rot_angle,'crop');   
                if wfs_transpose 
                    phiR=phiR';
                end       
               
                pad_dim     = (mdim-size(phiR))/2;
                phasek4     = padarray(phiR.*v_scale(aa),pad_dim,0,'both');
                a           = ifft2(pupil_amp.*exp(sqrt(-1).*phasek4));
                h           = abs(a).^2;
                hh          = hh + h;
            end  
            hh=hh./length(v_int);
            H(:,:,k)    = fft2(hh);
%             figure(8);imagesc(fftshift(hh));pause(0.5);
        end
        F       = sum(conj(H).*G,3)./(sum(abs(H).^2,3) + 1e-2);
        fval    = 0;
        for k=1:size(G,3)
            fval = fval + sum(sum(mask(:,:,k).*abs( G(:,:,k) - F.*H(:,:,k)).^2));
        end
        fx = real(ifft2(F));
        afx(:,:,aa,ss)  = fx;
%         v_max(aa,ss)    = max(fx(:));
         v_fit(aa,ss)    = fval;
%          figure(5);imagesc(phiR);pause(0.5);
%       %  figure(81);pcolor(v_size,v_scale,v_fit);colorbar;pause(0.2);shading interp
%         
%         figure(81);plot(v_scale(1:aa),v_fit(1:aa),'bsq-');colorbar;pause(0.2);
        
    end
%     ss
%     figure(11);   
%     set(gcf,'Color','w');
%     subplot(1,3,1);
% %     pcolor(v_size,v_scale,log10(v_fit));shading interp;colorbar
%     xlabel('Pupil diameter');ylabel('phase scaling');
%     set(gca,'FontSize',13);
%     title( [ 'wfs transp: ' num2str(wfs_transpose) '  angle: ' num2str(wfs_subap_rot_angle) ' > wind: ' num2str(vecs(1)) '  ' num2str(vecs(2))]);
% 
%     subplot(1,3,2);
%     pcolor(v_size,v_scale,v_max);shading interp
%     xlabel('Pupil diameter');ylabel('phase scaling');
%     set(gca,'FontSize',13);
%     title( [ 'wfs transp: ' num2str(wfs_transpose) '  angle: ' num2str(wfs_subap_rot_angle) ' > wind: ' num2str(vecs(1)) '  ' num2str(vecs(2))]);
%     colorbar
%     
%     pause(0.5);
end
% v_co = 129:129+255;
% 
% subplot(1,3,3);
% [aa,ss]=ind2sub(size(v_max),find(v_max==max(v_max(:))));[ v_scale(aa) v_size(ss) ]
% [aa,ss]=ind2sub(size(v_max),find(v_fit==min(v_fit(:))));[ v_scale(aa) v_size(ss) ]
% imagesc(afx(:,:,aa(1),ss(1)));colorbar;%title( [ 'wfs transpose: ' num2str(wfs_transpose) '  rotation: ' num2str(wfs_subap_rot_angle) ' ==> wind: ' num2str(vecs(1)) '  ' num2str(vecs(2))]);
% 


% 
% v_noise