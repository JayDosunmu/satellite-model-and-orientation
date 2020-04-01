function z = zernike(n,m,size)

% z = zernike(n,m,size)
%     Generates the (n,m) Zernike polynomial, where n>=m, on the inscribed circle
%     in a square raster of side 'size'. 
%
%     Inputs:
%	n	Radial order of polynomial; must be non-negative
%	m	Azimuthal order: in the range m = -n:2:n
%	size	Linear size of output raster
%
%     Output:
%	z	2D raster


ma = abs(m);
if (n<ma) | (n<0)
   error('n must be non-negative and n>=|m|');
end

s = (n-ma)/2;
if (s-floor(s)~=0)
   error('n-m must be even');
end

c1 = -((size-1)/2);
c2 = size+c1-1;

[x y] = meshgrid([c1:c2],[c1:c2]);
rho = sqrt(x.*x+y.*y)/size*2;
theta = atan2(y,x);

R = 0;
for k=0:s
   R = R+((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+ma)/2-k)*factorial(s-k)).*rho.^(n-2*k);
end

R = R.*circle(size,size,[(size+1)/2 (size+1)/2]);
if m<0
   z = sqrt(2*(n+1))*R.*sin(ma*theta);
elseif m>0
   z = sqrt(2*(n+1))*R.*cos(ma*theta);
else
   z = sqrt(n+1)*R;
end
