function y=smooth_moving(x,m) 

% y = smooth_moving(x,m)
%       x 	is the input vector (or matrix) to be smoothed. 
%       m 	is number of points to average over (best odd, but even works)
%       y 	is output vector of same length as x
%
% directly replaces the smooth.m from the Curvefitting toolbox
%

f=zeros(m,1)+1/m;
n=size(x,1);
isodd=bitand(m,1);
m2=floor(m/2);

y=filter(f,1,x);
y=y([zeros(1,m2-1+isodd)+m,m:n,zeros(1,m2)+n]);

for k=1:m2
   y(k) = mean(x(1:k*2-1));
   y(end-k+1) = mean(x(end-k*2+2:end));
end