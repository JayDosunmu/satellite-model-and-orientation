function zerns = zernikes(pdim,nz)
 
v=1:pdim;
icen = 1+(pdim-1)/2;
x = (v-icen);
x = x./abs(x(1));

[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
p = 0:(nz-1);
z = zeros(size(X));
y = zernfun2(p,r(idx),theta(idx));
zerns = zeros(pdim,pdim,nz);
for k=1:nz
    z(idx) = y(:,k);
    zerns(:,:,k) = z;
end


