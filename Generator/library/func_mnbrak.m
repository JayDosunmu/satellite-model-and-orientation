function [ax,bx,cx,fa,fb,fc] = func_mnbrak(func,ax,bx,x_in,d_in)

%1. Uses a user-defined func.m
%2. To bracket a minimum from initial guess of ax and bx; see p.393 in Numerical Recipes F77. 
%3. I added axt=x_in+ax.*d_in, etc, for use in multidimensional conjugate gradient 
%   calculation, where the search direction is d. In NR, this step is 
%   performed by the function f1dim on p.413.
%4. If you are following us, the outputs from these function are all scalars.
%5. Nothing too interesting here, unless you are hunting a bug.

GOLD=1.618034;
GLIMIT=100.;
TINY=1.e-20;

axt=x_in+ax.*d_in;
[fa]=feval(func,axt);
bxt=x_in+bx.*d_in;
[fb]=feval(func,bxt);

if fb > fa 
  dum=ax;
  ax=bx;
  bx=dum;
  dum=fb;
  fb=fa;
  fa=dum;
end

%first guess for c
cx=bx+GOLD.*(bx-ax);
cxt=x_in+cx.*d_in;
[fc]=feval(func,cxt);

while fb >= fc

     r=(bx-ax).*(fb-fc);
     q=(bx-cx).*(fb-fa);
     u=bx-((bx-cx).*q-(bx-ax).*r)./(2.*abs(max(abs(q-r),TINY)).*sign(q-r));
     
     if (2.*abs(max(abs(q-r),TINY)).*sign(q-r))==0
  %       3.134
     end
     ulim=bx+GLIMIT.*(cx-bx);
     if (bx-u).*(u-cx) > 0. 
        ut=x_in+u.*d_in;
        [fu]=feval(func,ut);
        if fu < fc  
           ax=bx;
           fa=fb;
           bx=u;
           fb=fu;
           break;
        elseif fu > fb 
           cx=u;
           fc=fu;
           break;
        end
        u=cx+GOLD.*(cx-bx);
        ut=x_in+u.*d_in;
        [fu]=feval(func,ut);
     elseif (cx-u).*(u-ulim) > 0. 
        ut=x_in+u.*d_in; 
        [fu]=feval(func,ut);
        if fu < fc 
           bx=cx;
           cx=u;
           u=cx+GOLD.*(cx-bx);
           fb=fc;
           fc=fu;
           ut=x_in+u.*d_in;
           [fu]=feval(func,ut);
        end
     elseif (u-ulim).*(ulim-cx) >= 0 
           u=ulim;
           ut=x_in+u.*d_in;
           [fu]=feval(func,ut);
     else
            u=cx+GOLD.*(cx-bx);
            ut=x_in+u.*d_in;
            [fu]=feval(func,ut);
     end
     ax=bx;
     bx=cx;
     cx=u;
     fa=fb;
     fb=fc;
     fc=fu;
end    