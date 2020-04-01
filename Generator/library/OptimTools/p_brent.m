function [xmin,fmin,iter,FunEvals] = p_brent(f, ax, bx, cx, x_in,d_in, tol)

%  Brents Method for Minimization.
%
%  [ xmin, fmin ] = brent(f, ax, bx, cx, x_in,d_in, tol)
%
%  Given a function f, and given a bracketing triplet of abscissas ax,
%  bx and cx, and fx(b) is less than both f(ax) and f(cx), this routine
%  isolates the minimum to a fractional precision of about TOL using
%  Brent's method. The abscissa of the minimum is returned in xmin, and
%  the minimum function value is returned in fmin.
%
%  References:
%
%  Brent, R.P. 1973, Algorithms for Minimization without Derivatives
%     (Englewood Cliffs, NJ: Prentice Hall), Chapter 5
%  Press, W.H. et al, 1992,  Numerical recipes in C: the art of scientific
%     computing (Cambridge University Press), Chapter 10.2
%
%  Description of variables:
%
%    a, b   the minimum is bracketed between a and b
%    x      point with the least function value found so far
%    w      point with the second least function value found so far
%    v      previous value of w
%    u      point at which the function was evaluated most recently
%    xm     midpoint between a and b (the function is not evaluated there)
%
%    Note: these points are not necessarily all distinct.
%
%    e      movement from best current value in the last iteration step
%    etemp  movement from the best current value in the second last
%           iteration step
%
%  General principles:
%
%    - The parabola is fitted trough the points x, v and w
%    - To be acceptable the parabolic step must
%        (i)  fall within the bounding interval (a,b)
%        (ii) imply a movement from the best current value x that is
%             *less* than half the movement of the *step before last*
%    - The code never evaluates the function less than a distance tol1
%      from a point already evaluated or from a known bracketing point
%
% Written by:
% --
% John L. Weatherwax                2004-12-11
%
% email: wax@alum.mit.edu
%
% Please send comments and especially bug reports to the
% above email address.
%
%-----


VERBOSE = 0;                % Print steps taken for the solution.
CGOLD   = (3-sqrt(5))/2;    % golden ratio
ITMAX   = 100;              % max number of iterations
ZEPS    = 1e-10;            % absolute error tolerance

% initialization
e = 0;
if ax < cx
 a = ax; b = cx;
else
 a = cx; b = ax;
end;
v = bx; w = v; x = w;

%x_in,d_in

fx = feval(f,x_in + x.*d_in); fv = fx; fw = fv;
FunEvals = 1;

for iter = 1:ITMAX,
 if( VERBOSE )
   fprintf(1, 'k=%4d, |a-b|=%e\n', iter, abs(a-b));
 end
 xm = 0.5*(a+b);
 tol1 = tol*abs(x) + ZEPS;
 tol2 = 2*tol1;
 % Stopping criterion: equivalent to: max(x-a, b-x) <= tol2
 if abs(x-xm) <= tol2-0.5*(b-a)
   xmin = x;
   fmin = fx;
   return
 end
 if abs(e) > tol1
   %
   % The second last move was sufficently large:
   % let's construct the parabolic fit
   %
   r = (x-w)*(fx-fv);
   q = (x-v)*(fx-fw);
   p = (x-v)*q - (x-w)*r;
   q = 2.0*(q-r);
   if q > 0, p = -p; end
   q = abs(q);
   etemp = e;
   e = d;
   if abs(p) >= abs(0.5*q*etemp) | p <= q*(a-x) | p >= q*(b-x)
     %
     % The parabolic fit did not meet the reqirements above:
     % use a golden section refinement instead.
     %
     % Print an explanation of what condition failed:
     %
     if( VERBOSE )
       if abs(p) >= abs(0.5*q*etemp)
        % disp('abs(p) >= abs(0.5*q*etemp)')
       end
       if p <= q*(a-x)
         %disp('p <= q*(a-x)')
       end
       if p >= q*(b-x)
         %disp('p >= q*(b-x)')
       end
     end
     %
     % Set e in such a way, that a parabolic fit is possible in the
     % next iteration.
     %
     if x >= xm
       e = a-x;
     else
       e = b-x;
     end
     d = CGOLD*e;
     if( VERBOSE )
       disp('golden section step');
     end
   else
     %
     % the parabolic fit meets the above requirements: use iteration
     %
     d = p/q;
     u = x+d;
     if u-a < tol2 | b-u < tol2
       %
       % the parabolic fit is too close to the boundaries of the
       % bracketing interval: put new point tol1 apart from the
       % midpoint
       %
       d = tol1*sign(xm-x);
       if( VERBOSE ) disp('too close to interval boundaries'); end
     else
       if( VERBOSE ) disp('parabolic step'); end
     end
   end
 else
   %
   % The second last move was not sufficently large: use golden section
   % set e in such a way, that a parabolic fit is possible in the next
   % iteration.
   %
   if x >= xm
     e = a-x;
   else
     e = b-x;
   end
   d = CGOLD*e;
   if( VERBOSE )
     disp('abs(e) > tol1');
     disp('golden section step');
   end
 end
 if abs(d) >= tol1
   u = x+d;
 else
   %
   % u is too close to x: put u farer apart
   %
   u = x + tol1*sign(d);
   if( VERBOSE ) disp('u is too close to x'); end;
 end
 fu = feval(f,x_in + u.*d_in);    % finally evaluate the function
 FunEvals = FunEvals + 1;
 if fu <= fx
   %
   % u is better than x: x becomes new boundary
   %
   if u >= x
     a = x;
   else
     b = x;
   end
   v = w; w = x; x = u;
   fv = fw; fw = fx; fx = fu;
 else
   %
   % x is better than u: u becomes new boundary
   %
   if u < x
     a = u;
   else
     b = u;
   end
   if fu <= fw | w == x
     v = w; w = u;
     fv = fw; fw = fu;
   elseif fu <= fv | v == x | v == w
     v = u;
     fv = fu;
   end
 end
end
error('too many iterations in brent');