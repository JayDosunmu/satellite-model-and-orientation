function [ phasep, zc ] = sub_modes(phase_maps,Zmodes,Zfrac,z_basis,mask)

% subtract of zernike modes from a phase screen
% 
% [ phasep, coeff ] = sub_modes( phases, Zmodes, Zfrac, Zbasis )
%
if nargin == 4
    mask = ones(size(phase_maps,1),size(phase_maps,2));
end

phasep      = zeros(size(phase_maps));
phi_energy  = zeros(Zmodes,size(phase_maps,3));
zc          = zeros(Zmodes,size(phase_maps,3));

for k=1:size(phase_maps,3)
    tmp     = phase_maps(:,:,k);
    for q=1:Zmodes
        c_q      = sum(sum(mask.*(tmp.*z_basis(:,:,q))));
        c_area   = sum(sum(mask.*(z_basis(:,:,q)).^2));
        a(q)        = c_q/c_area;  
        zc(q,k)     = a(q);
         amp         = a(q);% c_q/c_area; 
      
        if (q <= Zmodes)
            tmp     = tmp - Zfrac*amp.*(mask.*z_basis(:,:,q));   
        end      
    %    tmp=tmp.*pupil_mask;
    end  
    phasep(:,:,k) = tmp.*z_basis(:,:,1);  
end
