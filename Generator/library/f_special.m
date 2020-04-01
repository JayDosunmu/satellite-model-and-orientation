function h = f_special(ctype,sz,sigx)


[x,y]   = meshgrid(-(sz/2):1:(sz/2)-1,-(sz/2):1:(sz/2)-1);
rad2    = x.^2+y.^2;
h       = exp(-rad2./(2*sigx^2))./(2*pi*sigx^2);
