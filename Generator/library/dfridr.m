function [ err, ans, a ] = dfridr( func_name, x, pp, hh) 

CON = 1.4;
CON2 = CON*CON;
SAFE = 2.0;

global cost_params

xp1=x;
xp2=x;
a=zeros(10,10);
x1=x;x2=x;
x1(pp)=x1(pp)+hh;
x2(pp)=x2(pp)-hh;
a(1,1)=(feval(func_name,x1) - feval(func_name,x2))/(2*hh);
err=1e30;
for i=2:20
    hh = hh/CON;
    xp1=x;
    xp1(pp) = xp1(pp) - hh;
    xp2(pp) = xp2(pp) + hh;
    a(1,i)=(feval(func_name,xp2) - feval(func_name,xp1))/(2*hh);
    fac = CON2;
    for j=2:i
        a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1.0);
        fac=CON2*fac;
        errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)));
        if ( errt <= err )
            err=errt;
            ans=a(j,i);
        end
    end
    if ( abs(a(i,i)-a(i-1,i-1)) >= SAFE*err )
        break
    end
    xp1(pp) = xp1(pp) + hh;
    xp2(pp) = xp2(pp) - hh;    
end