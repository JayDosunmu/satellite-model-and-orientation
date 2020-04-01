function [ p, IterInfo] = p_NLCG(p, OPTIM_params, fval0)
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
%   objectç - current estimate of the object. 
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
%                   FunVals - vector containg the function values at each
%                             iteration.
%                  FunEvals - vector containing the number of function 
%                             evaluations needed for each iteration 
%                 GradEvals - vector containing the number of gradient 
%                             evaluations needed for each iteration 
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

global cost_params
% object = reshape(p(cost_params.idx_obj_vars).^2,[cost_params.mdim cost_params.mdim]);

%
% First extract some information from the OPTIM_params structure.
%
func_name           = OPTIM_params.func_name;
func_and_dfunc_name = OPTIM_params.func_and_dfunc_name;
MaxIters            = OPTIM_params.NLCG_MaxIters;
RelFunChangeTol     = OPTIM_params.NLCG_RelFunChangeTol;
RelGradTol          = OPTIM_params.NLCG_RelGradTol;
IterPrint           = OPTIM_params.IterPrint;


%
% If IterPrint = 0, then we should print information at each iteration
%
if IterPrint == 0
  echo_out(' Iter  FuncEvals GradEvals    F(X)/N        |F(X)-F(Xold)|/|F(Xold)|   ||G(X)||/N           Time (seconds)     \n');
  echo_out('------ --------- --------- ---------------- ------------------------  ----------------     --------------- \n');
end
tic;

%
%  To call the func and dfunc, it will be easier to organize the precise
%  inputs the function needs, and just pass this one structure.
%  Here we put together this input structure.
%
% fun_input.DORA_input_data = DORA_input_data;
% fun_input.WFS_params      = WFS_params;
% fun_input.PSF_params      = PSF_params;
% fun_input.FFM_options     = FFM_options;
% fun_input.OPTIM_params    = OPTIM_params;
% fun_input.OBJ_params      = OBJ_params;
% fun_input.SparseMatrices  = SparseMatrices;

%
%  We will use a function handle to evaluate the function and gradient,
%  so here we define those function handles.
%
func = eval(['@(p)',func_name,'(p)']);
func_and_dfunc  = eval(['@(p)',func_and_dfunc_name,'(p)']);

%
%  We need a tolerance for finding a line search using Brent's method
%  (see below).  It might be worth investigating what should be the
%  best choice for this parameter.
%
Brent_TOL = 2e-4;
 
%
%  We need to compute an initial gradient, and the initial direction
%  d = negative gradient.
%
%g = dfunc(p);
[fval, g] = feval(func_and_dfunc, p);
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
    %[ a_hat, H_hat, aAu, aPu ]      = feval(cost_params.PSF_func,p);
   % object = reshape(p(cost_params.idx_obj_vars).^2,[cost_params.mdim cost_params.mdim]);
  %  GH = fft2(object).*H_hat{cost_params.curr_chann}.H_hat(:,:,1);
 %   echo_out(sprintf('%6d %9d %9d %16.8e   %16.8e       %16.8e      %8.2f   %10.5e   %10.5e   \n', 0, FunEvals(1), GradEvals(1), fval, RelFunChange(1), RelGradNorms(1), IterTimes(1),sum(object(:)),GH(1,1)./sum(sum(cost_params.gg2{cost_params.curr_chann}(:,:,1)))));
  echo_out(sprintf('%6d %9d %9d %16.8e   %16.8e       %16.8e      %8.2f    \n', 0, FunEvals(1), GradEvals(1), fval, RelFunChange(1), RelGradNorms(1), IterTimes(1)));
end
tic;

if ~isempty(OPTIM_params.out_func)
    feval(OPTIM_params.out_func,p,0,OPTIM_params);
end

for ii = 1:MaxIters
    cost_params.ii = cost_params.ii +1;
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
  if ii == 1
    ax = 0;
    bx = 1;
  end
  [ax,bx,cx,fa,fb,fc,BRAK_FunEvals] = p_func_mnbrak(func,ax,bx,p,d);
  FunEvalsii = FunEvalsii + BRAK_FunEvals;
  %
  % There are built-in functions in MATLAB called "lower" and "upper", so
  % don't use those as variable names here.
  %
  if (ax > cx)
    lower_bkt=cx;
    upper_bkt=ax;
  else
    lower_bkt=ax;
    upper_bkt=cx;
  end               
  [alpha, ~, ~, BRENT_FunEvals] = p_brent(func, lower_bkt, bx, upper_bkt,p,d, Brent_TOL);
  FunEvalsii = FunEvalsii + BRENT_FunEvals;
  ax = 0;
  bx = alpha;
  
  %
  %  Note that I think the above function brent.m is essentially
  %  the same as what MATLAB's built-in fminbnd function does.
  %  It would be interesting to compare this.
  %  It might also be worth changing this exact line search to an
  %  inexact line search based on the strong Wolfe condition.
  %
  %p    = p + alpha*xi;
  p = p + alpha*d;

  %
  % Save the previous function value and gradint, and update:
  %
  fvalp = fval;
  g_old = g;
  [fval, g] = feval(func_and_dfunc, p);
  FunEvalsii = FunEvalsii + 1;
  GradEvalsii = GradEvalsii + 1;
           
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
%   [ a_hat, H_hat, aAu, aPu ]      = feval(cost_params.PSF_func,p);
%   object = reshape(p(cost_params.idx_obj_vars).^2,[cost_params.mdim cost_params.mdim]);
%     GH = fft2(object).*H_hat{cost_params.curr_chann}.H_hat(:,:,1);
  if IterPrint == 0
    echo_out(sprintf('%6d %9d %9d %16.8e   %16.8e       %16.8e      %8.2f    \n', ii, FunEvals(ii+1), GradEvals(ii+1), fval, RelFunChange(ii+1), RelGradNorms(ii+1), IterTimes(ii+1)));
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
  if ~isempty(OPTIM_params.out_func)
  	feval(OPTIM_params.out_func,p,ii,OPTIM_params);
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
    echo_out('^^^^^^\n');
    echo_out('NLCG STOP\n');
    echo_out('MaxIter\n');
    echo_out(sprintf('%6d\n',MaxIters));
  elseif StopFlag == 2 
    echo_out('                                               ^^^^^^^^^^^^^^^^\n');
    echo_out('                                                  NLCG STOP\n');
    echo_out('                                                RelFunChangeTol\n');
    echo_out(sprintf('                                              %16.8e\n',RelFunChangeTol));
  else
    echo_out('                                                                       ^^^^^^^^^^^^^^^^\n');
    echo_out('                                                                          NLCG STOP\n');
    echo_out('                                                                        RelGradNormTol\n');
    echo_out(sprintf('                                                                     %16.8e\n',RelGradTol));
  end
end


