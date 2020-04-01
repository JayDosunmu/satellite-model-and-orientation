function y = p_Sinc(x)
%
% This simple function just evaluates 
%
%   y = sinc(x) = / sin(x)/x, if x ~= 0
%                 \ 1,        if x == 0
%
% at an array of x-values.  
%
y = ones(size(x));
y(x~=0) = sin(x(x~=0))./x(x~=0);
