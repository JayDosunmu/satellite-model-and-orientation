function mask = make_circle_mask(sz,rad)

% mask = dora_make_mask(sz,rad)

c1 = -((sz-1)/2);
c2 = sz+c1-1;

[x,y]   = meshgrid([c1:c2],[c1:c2]);
rho     = sqrt(x.*x+y.*y);



theta = atan2(y,x);

mask=circle(sz,sz,[(sz+1)/2 (sz+1)/2]);
mask(rho>rad)=0.0;



