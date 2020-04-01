function h = func_gaussian(mdim,sig2)

if sig2>0
    mdim2   = mdim/2;
    [x,y]   = meshgrid(-mdim2:1:(mdim2-1),-mdim2:1:(mdim2-1));rad2=x.^2 + y.^2;       
    h       = exp(-rad2./(2*sig2^2))./(2*pi.*sig2^2);
else
    h = zeros(mdim);
    h(mdim/2+1,mdim/2+1)=1.0;
end

