function dpp = fnder_spline(pp, d)
%
%  Find the d'th derivative of a cubic spline.
%  If d is not specified, find the first derivative.
%

if nargin == 1
  d = 1;
end

if d < 0
  error('d must be >= 0')
end
switch d
  case 0
    dpp = pp
  case 1
    [xb, pp_coef] = unmkpp(pp);
    dpp_coef = [3*pp_coef(:,1), 2*pp_coef(:,2), pp_coef(:,3)];
    dpp = mkpp(xb, dpp_coef);
  case 2
    [xb, pp_coef] = unmkpp(pp);
    dpp_coef = [6*pp_coef(:,1), 2*pp_coef(:,2)];
    dpp = mkpp(xb, dpp_coef);
  case 3
    [xb, pp_coef] = unmkpp(pp);
    dpp_coef = 6*pp_coef(:,1);
    dpp = mkpp(xb, dpp_coef);
  otherwise
    [xb, pp_coef] = unmkpp(pp);
    dpp_coef = zeros(size(pp_coef,1));
    dpp = mkpp(xb, dpp_coef);
end


