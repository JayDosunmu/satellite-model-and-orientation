function [ new_pupil_mask, new_pupil_amp, gg2, dm2 ] = upsample_data2newgrid(native_grid_sz, PupilDiam1,PupilDiam2, gg, dm, br_factor  )

% [ new_pupil_mask, new_pupil_amp, gg2, dm2 ] = upsample_data2newgrid( native_grid_sz, PupilDiam1,PupilDiam2, g, dm, br_factor  )
%
%
grid_sz = native_grid_sz*br_factor;

newPupil_Diam1 = br_factor*PupilDiam1;
newPupil_Diam2 = br_factor*PupilDiam2;
c1=circle(grid_sz,newPupil_Diam1);
c2=circle(grid_sz,newPupil_Diam2);
new_pupil_mask = c1-c2;

g_sz0       = size(gg,1);
g_pad_sz    = (br_factor.*g_sz0  - g_sz0 )/2;
gg2          = padarray(gg,[g_pad_sz g_pad_sz],0,'both');

dm2         = padarray(dm,[g_pad_sz g_pad_sz],0,'both');


new_pupil_amp  = new_pupil_mask./max(new_pupil_mask(:));
new_pupil_amp   = ( grid_sz.*new_pupil_amp )./sqrt(sum(sum(new_pupil_amp.^2)));
a = ifft2(new_pupil_amp); h = abs(a).^2;
fprintf('INFO: PSF computed on grid size %4.0f x %4.0f has volume %4.0f \n',size(new_pupil_amp,1),size(new_pupil_amp,2),sum(h(:)));

 
