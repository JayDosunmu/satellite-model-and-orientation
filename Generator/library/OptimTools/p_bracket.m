function [a, b, FunEvals] = p_bracket(FUN, c, x, d, phi0, phi_der0, LineSearch_params, FunEvals)
%
%

epsk = LineSearch_params.epsk;
rho  = LineSearch_params.rho;
MaxSteps = 4; %LineSearch_params.bracketMaxSteps;

c_vec = zeros(MaxSteps, 1);
phic_vec = zeros(MaxSteps,1);
phic_vec(1) = phi_der0;

c_vec(2) = c;

for j = 2:MaxSteps
  xc = x + c_vec(j)*d;
  [~, gc] = feval(FUN, xc);
  FunEvals = FunEvals + 1;
  phi_derc = d'*gc;
  if phi_derc >= 0
    b = c_vec(j);
    idx = find(phic_vec(1:j-1) <= phi0+epsk, 1, 'last');
    a = c_vec(idx);
    return
  elseif (phi_derc < 0  && phi_derc > phi0+epsk)
    [a, b, FunEvals] = p_update0(FUN, a, b, phi0, phi_der0, x, d, LineSearch_params, FunEvals);
    return
  else
    c(j+1) = rho*c(j);
  end
end
if j == MaxSteps
  error('bracket did max steps')
end
