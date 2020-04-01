function c = p_initial0(x, fval, g, LineSearch_params)
%
%  This should be used only for first iteration
%
%  x    = initial guess at unknowns
%  fval = initial function value
%  g    = initial gradient 
%  LineSearc_params = structure containing stuff needed for the line
%                     search. 

psi0 = LineSearch_params.psi0;

%
%  Need to check if some values are approximately zero. So use
%  MATLAB's eps for this (double precision is approximately 2.2204e-16).
%
zeroTol = eps;

if norm(x) > zeroTol
  c = psi0*max(abs(x))/max(abs(g));
elseif abs(fval) > zeroTol
  c = psi0*abs(fval)/(g'*g);
else
  c = 1;
end
