function [x_c, y_c, y] = FindCorner(x, y)
%
%  Find the first point of maximum curvature of a set of data.
%

x = x(:);
y = y(:);

% First performa a local smoothing of y-data
y = smooth_moving(y,4);
%y = smooth(y,4);

%
% Now fit a order 4 spline to the smoothed data
%
n = length(x);
pp = spline(x,y);

%
%  Extract abscissa and ordinate splines and differentiat them.
%  Compute the function values as default in spleval.
P     = spleval(pp);  dpp   = fnder_spline(pp);
D     = spleval(dpp); ddpp  = fnder_spline(pp,2);
DD    = spleval(ddpp);
ppx   = P(1,:);       ppy   = P(2,:);
dppx  = D(1,:);       dppy  = D(2,:);
ddppx = DD(1,:);      ddppy = DD(2,:);

% Compute the corner of the discretized .spline curve via max. curvature.
% Define curvature = 0 where both dppx and dppy are zero.
k1    = dppx.*ddppy - ddppx.*dppy;
k2    = (dppx.^2 + dppy.^2).^(1.5);
I_nz  = find(k2 ~= 0);
kappa = zeros(1,length(dppx));
kappa(I_nz) = -k1(I_nz)./k2(I_nz);
[kmax,ikmax] = max(kappa);
x_corner = ppx(ikmax); 
y_corner = ppy(ikmax);

% Locate the point on the discrete L-curve which is closest to the
% corner of the spline curve. If the curvature is negative everywhere,
% then define the leftmost point of the L-curve as the corner.
if (kmax < 0)
  x_c = x(n); y_c = y(n);
else
  index = find(x < x_corner & y < y_corner);
  if (length(index) > 0)
    [dummy,rpi] = min((x(index)-x_corner).^2 + (y(index)-y_corner).^2);
    rpi = index(rpi);
  else
    [dummy,rpi] = min((x-x_corner).^2 + (y-y_corner).^2);
  end
  x_c = x(rpi); y_c = y(rpi);
end
% 
% figure(1), clf
% plot(x, y, 'LineWidth', 2)
% hold on
% plot(x_c, y_c, 'ro', 'LineWidth', 2)