function check_derivatives( func_name, dfunc_name, p, vdx, max_checks)

global cost_params

if nargin == 4
    max_checks = 4;
end

% fprintf('Derivative check \n');
% fprintf('-------------------------------------------------------------\n');
% fprintf(['func_name  : ' func_name '\n' ]);
% fprintf(['dfunc_name : ' dfunc_name '\n\n' ]);
 
[ fval, grad ]       = feval(dfunc_name,p); 

delta_x = [ 4 3 2 1 0.1 0.01];

fprintf('initial step size used\n');
fprintf('%8.7g ',delta_x');

fprintf('\ncomparing analytical and numerical derivatives ... \n\n');
fprintf('                                                \t\t\t\t\t\t\t lowest error\n');        
fprintf('pixel  \t var value \t\t step \t error  \t numerical \t rel. error  \t analytic \t ratio\n');

% cost_params.Hupdate=1;
% cost_params.Fupdate=1;
% 

v_pixels            = vdx; %(1:2:end);

checks=0;
keep_going=1;
while checks < max_checks

       % pp=max(1,floor(2.8e4.*rand));

    pp=max(1,floor(length(v_pixels).*rand));
    if abs(grad(v_pixels(pp))) > 1e-15
        fprintf('%5.0f\t%8.5g\t',v_pixels(pp),p(v_pixels(pp)));
        zerr=zeros(length(delta_x),1);zfd=zerr;
        for zp=1:length(delta_x)
            [ err, fd, a ] = dfridr(func_name,p,pp,delta_x(zp));
            zerr(zp)=err;zfd(zp) = fd;
        end
        idx=find(zerr==min(zerr(:))); [ err, fd, a ] = dfridr(func_name,p,v_pixels(pp),delta_x(idx(1)));
        fprintf('%8.5g \t',delta_x(idx(1)));
        fprintf('%8.5g \t %8.5g\t',err, fd);
        fprintf('%8.5g\t %8.5g\t%8.5g\n',abs(err/fd),grad(v_pixels(pp)),grad(v_pixels(pp))/fd);
        %    end
        checks = checks + 1;
        
    end
end
