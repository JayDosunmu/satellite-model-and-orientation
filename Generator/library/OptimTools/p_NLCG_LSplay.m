function [ p, IterInfo] = p_NLCG_LSplay(p, object, DORA_input_data, WFS_params, PSF_params, OBJ_params, OPTIM_params, FFM_options, SparseMatrices, fval0)
%
% Nonlinear conjugate gradient (NLCG) method, used to update the PSF 
% parameters, which are the high resolution composite phase gradients on 
% each atmospheric layer.
%
% This implementation is a translation from Numerical Recipes.  It miht
% be worth changing this, for example to use a Strong Wolfe condition for
% the line search, instead of the current exact line search we are doing.
%
% Input:  
%   p      - initial guess of the PSF paramters -- these are the
%            unknowns over which we want to optimize.
%   object - current estimate of the object.
%   Structures:
%     The folowing strucutres contain a variety of parameters, which are
%     needed for the function evaluations.  See other codes for the 
%     description of these parameters, as well as 
%     p_GradOptFun_and_DFun_poly.m for more information.
%          - DORA_input_data
%          - WFS_params
%          - PSF_params
%          - OBJ_params
%          - OPTIM_params
%          - FFM_options
%   SparseMatrices - cell array containing the sparse matrices from the
%                    FFM reconstruction.  These are needed in the
%                    function and gradient evaluations.
%
%   Optional Input:
%    fval0 - initial function value -- used to initialize RelFunChange.
%            Default to is to initialize RelFunChange(1) = 1;
%
%  Output:
%    p - update of the PSF parameters
%    IterInfo - structure that contains the following information about
%               the iterations:
%                  StopFlag - indicates why the iteration terminated:
%                             StopFlag = 0 
%                               ==> MaxIters reached
%                             StopFlag = 1 
%                               ==> norm(gradient)/length(p) is less than 
%                                   OPTIM_params.NLCG_RelGradTol
%                             StopFlag = 2 
%                               ==> |F(p)-F(p_old)|/|F(p_old)| is less than
%                                   OPTIM_params.NLCG_RelFunChangeTol
%                  FunEvals - vector containing the number of iterations
%                             needed for each iteration (time for 
%                             initialization is the first entry, time for
%                             iteration 1 is the second entry, etc.)
%                  FunEvals - vector containing the number of iterations
%                             needed for each iteration (time for 
%                             initialization is the first entry, time for
%                             iteration 1 is the second entry, etc.)
%              RelGradNorms - vector containing norm(gradient)/length(p)
%                             for each iteration (first entry is 1, to 
%                             correspond to initialization step, second
%                             entry is for iteration 1, etc.)
%              RelFunChange - vector containing |F(p)-F(p_old)|/|F(p_old)|
%                             for each iteration (first entry is 1, to 
%                             correspond to initialization step, second
%                             entry is for iteration 1, etc.)
%                 IterTimes - vector containing time needed for each
%                             iteration (first entry is for the 
%                             initialization step, second entry is for 
%                             iteration 1, etc.)
%                 

%LineSearch_params.delta = 0.1;
%LineSearch_params.sigma = 0.9;
%LineSearch_params.epsk  = 1e-6;
%LineSearch_params.omega = 1e-3;
%LineSearch_params.theta = 0.5;
%LineSearch_params.gamma = 0.66;
%LineSearch_params.rho   = 5;
LineSearch_params.ftol   = 1e-4;
LineSearch_params.gtol   = 1e-2;
LineSearch_params.xtol   = 1e-4;
LineSearch_params.stpmin = 1e-15;
LineSearch_params.stpmax = 1e15;
LineSearch_params.maxfev = 20;
LineSearch_params.psi0   = 1e-2;
LineSearch_params.psi1   = 1e-1;
LineSearch_params.psi2   = 2;

%
% First extract some information from the OPTIM_params structure.
%
%func_name       = OPTIM_params.func_name;
%dfunc_name      = OPTIM_params.dfunc_name;
func_and_dfunc_name = OPTIM_params.func_and_dfunc_name;
MaxIters        = OPTIM_params.NLCG_MaxIters;
RelFunChangeTol = OPTIM_params.NLCG_RelFunChangeTol;
RelGradTol      = OPTIM_params.NLCG_RelGradTol;
IterPrint       = OPTIM_params.IterPrint;
LogFile_fid     = OPTIM_params.LogFile_fid;

%
% If IterPrint = 0, then we should print information at each iteration
%
if IterPrint == 0
  fprintf(LogFile_fid,' Iter  FuncEvals GradEvals    F(X)/N        |F(X)-F(Xold)|/|F(Xold)|   ||G(X)||/N           Time (seconds)     \n');
  fprintf(LogFile_fid,'------ --------- --------- ---------------- ------------------------  ----------------     --------------- \n');
  tic;
end

%
%  To call the func and dfunc, it will be easier to organize the precise
%  inputs the function needs, and just pass this one structure.
%  Here we put together this input structure.
%
fun_input.DORA_input_data = DORA_input_data;
fun_input.WFS_params      = WFS_params;
fun_input.PSF_params      = PSF_params;
fun_input.FFM_options     = FFM_options;
fun_input.OPTIM_params    = OPTIM_params;
fun_input.OBJ_params      = OBJ_params;
fun_input.SparseMatrices  = SparseMatrices;

%
%  We will use a function handle to evaluate the function and gradient,
%  so here we define those function handles.
%
%func  = eval(['@(p)',func_name,'(p, object, fun_input)']);
%dfunc = eval(['@(p)',dfunc_name,'(p, object, fun_input)']);
func_and_dfunc  = eval(['@(p)',func_and_dfunc_name,'(p, object, fun_input)']);

%
%  We need a tolerance for finding a line search using Brent's method
%  (see below).  It might be worth investigating what should be the
%  best choice for this parameter.
%
%Brent_TOL = 2e-4;
 
%
%  We need to compute an initial gradient, and the initial direction
%  d = negative gradient.
%
[fval, g] = feval(func_and_dfunc, p);
%g = dfunc(p);
d = -g;

%  And we need to compute an initial function evaluation.
%fval = func(p);
%fval0 = fval;

%
%  Before starting the iterations, we need to initialize several things:
%
n_params = length(p);
StopFlag = 0;
IterTimes = zeros(MaxIters+1,1);
FunEvals = zeros(MaxIters+1,1);
FunVals = zeros(MaxIters+1,1);
GradEvals = zeros(MaxIters+1,1);
RelGradNorms = zeros(MaxIters+1,1);
RelFunChange = zeros(MaxIters+1,1);
FunVals(1) = fval;
FunEvals(1) = 1;
GradEvals(1) = 1;
FunEvalsii = 0;
GradEvalsii = 0;
RelGradNorms(1) = norm(g)/n_params;
if isempty(fval0)
  RelFunChange(1) = 1;
else
  RelFunChange(1) = abs(fval-fval0)/abs(fval0);
end
IterTimes(1) = toc;
%
%  If IterPrint = 0, print to the screen information for  the initialization
%  step (iteration 0).
%
if IterPrint == 0
  fprintf(LogFile_fid,'%6d %9d %9d %16.8e   %16.8e       %16.8e      %8.2f \n', 0, FunEvals(1), GradEvals(1), fval/n_params, RelFunChange(1), RelGradNorms(1), IterTimes(1));
end
tic;

for ii = 1:MaxIters
  %
  %  In nonlinear conjugate gradient, the first thing that needs to
  %  be done is to use a line search method to compute a step
  %  step length.  Basically, the idea is to do the following:
  %  find alpha to minimize:
  %      f(alpha) = func_name(p + alpha*d)
  %  where p is the parameter vector and d is the current step
  %  direction.  Note that, although func_name is a multivariable
  %  function, for the line search, we only consider f(alpha) as
  %  a function of a single variable.
  %
  %  There are many ways to find this minimum ... here we are using
  %  "Brent's method", which needs a bracket around where the 
  %  minimum is located.  
  %  
  %alpha = 1;
  %if ii == 1;
    alpha = p_initial0(p, fval, g, LineSearch_params);
  %else
  %  phi0 = fval;
  %  phi_der0 = d'*g;
  %  [alpha, FunEvals_initial] = p_initial(func_and_dfunc, p, d, phi0, phi_der0, alpha, LineSearch_params);
  %  FunEvalsii = FunEvalsii + FunEvals_initial;
  %  GradEvalsii = GradEvalsii + FunEvals_initial;
  %end
  fvalp = fval;
  g_old = g;
  [p,fval,g,alpha,info,LS_FunEvals] = p_cvsrch(func_and_dfunc,p,fval,g,d,alpha,LineSearch_params);
  FunEvalsii = FunEvalsii + LS_FunEvals;
  GradEvalsii = GradEvalsii + LS_FunEvals;
  %[alpha_new, FunEvalsiii]
  %if ii == 1
  %  ax = 0;
  %  bx = 1;
  %end
  %[ax,bx,cx,fa,fb,fc,FunEvalsii] = p_func_mnbrak(func,ax,bx,p,d,FunEvalsii);
  %
  % There are built-in functions in MATLAB called "lower" and "upper", so
  % don't use those as variable names here.
  %
  %if (ax > cx)
  %  lower_bkt=cx;
  %  upper_bkt=ax;
  %else
  %  lower_bkt=ax;
  %  upper_bkt=cx;
  %end               
  %[alpha, ~, ~, FunEvalsii] = p_brent(func, lower_bkt, bx, upper_bkt,p,d, Brent_TOL, FunEvalsii);
  %ax = 0;
  %bx = alpha
  
  %[alpha, alpha_new]
  %[FunEvalsii, FunEvalsiii]
  %
  %  Note that I think the above function brent.m is essentially
  %  the same as what MATLAB's built-in fminbnd function does.
  %  It would be interesting to compare this.
  %  It might also be worth changing this exact line search to an
  %  inexact line search based on the strong Wolfe condition.
  %
  %p    = p + alpha*d;

  %
  % Save the previous function value, and compute the current
  % function value:
  %
  %fvalp = fval;
  %fval = func(p);
  %FunEvalsii = FunEvalsii + 1;

  % udate the gradient
  %g_old = g;
  %g = dfunc(p);
  %FunEvalsii = FunEvalsii + 1;
           
  % Use Polak-Ribiere CG 
  gg = sum(g_old(:).*g_old(:));
  dgg = sum((g(:)-g_old(:)).*g(:));
  beta    = dgg/gg;
  %  If you want to use Fletcher-Reeves CG, then uncomment the 
  %  next few lines, and comment out the Polak-Ribiere lines.
  %
  % gg      = sum(g_old(:).*g_old(:));
  % dgg     = sum(g(:).*g(:));
  % beta    = dgg/gg;
  %
  d = -g + beta*d; 
  
  %
  % Update iteration information, and check for convergence
  %
  RelGradNorms(ii+1) = norm(g)/n_params;
  RelFunChange(ii+1) = abs(fval-fvalp)/abs(fvalp);
  FunVals(ii+1) = fval;
  FunEvals(ii+1) = FunEvalsii;
  GradEvals(ii+1) = GradEvalsii;
  FunEvalsii = 0;
  GradEvalsii = 0;
  IterTimes(ii+1) = toc;
  %
  %  If IterPrint = 0, print the iteration information to the screen
  %
  if IterPrint == 0
    fprintf(LogFile_fid,'%6d %9d %9d %16.8e   %16.8e       %16.8e      %8.2f \n', ii, FunEvals(ii+1), GradEvals(ii+1), fval/n_params, RelFunChange(ii+1), RelGradNorms(ii+1), IterTimes(ii+1));
  end
  tic;
  if RelGradNorms(ii+1) <= RelGradTol
    StopFlag = 1;
    break
  end
  if RelFunChange(ii+1) <= RelFunChangeTol
    StopFlag = 2;
    break
  end
end
IterInfo.StopFlag     = StopFlag;
IterInfo.FunEvals     = FunEvals(1:ii+1);
IterInfo.GradEvals    = GradEvals(1:ii+1);
IterInfo.FunVals      = FunVals(1:ii+1);
IterInfo.RelGradNorms = RelGradNorms(1:ii+1);
IterInfo.RelFunChange = RelFunChange(1:ii+1);
IterInfo.IterTimes    = IterTimes(1:ii+1);
if IterPrint == 0
  if StopFlag == 0
    fprintf(LogFile_fid,'^^^^^^\n');
    fprintf(LogFile_fid,'NLCG STOP\n');
    fprintf(LogFile_fid,'MaxIter\n');
    fprintf(LogFile_fid,'%6d\n',MaxIters);
  elseif StopFlag == 2 
    fprintf(LogFile_fid,'                                               ^^^^^^^^^^^^^^^^\n');
    fprintf(LogFile_fid,'                                                  NLCG STOP\n');
    fprintf(LogFile_fid,'                                                RelFunChangeTol\n');
    fprintf(LogFile_fid,'                                              %16.8e\n',RelFunChangeTol);
  else
    fprintf(LogFile_fid,'                                                                      ^^^^^^^^^^^^^^^^\n');
    fprintf(LogFile_fid,'                                                                         NLCG STOP\n');
    fprintf(LogFile_fid,'                                                                       RelGradNormTol\n');
    fprintf(LogFile_fid,'                                                                     %16.8e\n',RelGradTol);
  end
end


