function [phase] = genphase(num_increases)
%  Generates a phase screen with Kolmogorov statistics.
%
%  function [phase] = genphase(num_increases)
%
%  Inputs       num_increases   - the number of times the original 15x15
%                                 screen is to be interpolated
%
%  Outputs      phase           - the interpolated phase screen with the
%                                 edges removed
%
%  Rachel Johnston and Cressida Harding
%  Department of Electrical and Electronic Engineering
%  University of Canterbury
%  Christchurch
%  New Zealand
%
%  Code may be copied and used freely provided acknowledgement of the original source
%  is maintained.
%
% First load the singular vectors and values of a 15x 15 turbulence
% [u,d,v]. Code should be modified to read the file outside the function
% and pass u,d,v as parameters if greater speed is required.
%
load svectors4
%
residual = 0.08441735664383*3^(5/3);

inter1 = [...
-0.00166566566045,  0,-0.03407545836221, 0,-0.03407545836221, 0,-0.00166566566045;...
               0,  0,                0, 0,                0, 0,                0;...
-0.03407545836221,  0, 0.31981657662220, 0, 0.31981657662220, 0,-0.03407545836221;...
               0,  0,                0, 1,                0, 0,                0;...
-0.03407545836221,  0, 0.31981657662220, 0, 0.31981657662220, 0,-0.03407545836221;...
               0,  0,                0, 0,                0, 0,                0;...
-0.00166566566045,  0,-0.03407545836221, 0,-0.03407545836221, 0,-0.00166566566045];

inter2 = [ ...
               0,                0,                0,-0.00166566566045,                0,                0,                0;...
               0,                0,-0.03407545836221,                0,-0.03407545836221,                0,                0;...
               0,-0.03407545836221,                0, 0.31981657662220,                0,-0.03407545836221,                0;...
-0.00166566566045,                0, 0.31981657662220,                1, 0.31981657662220,                0,-0.00166566566045;...
               0,-0.03407545836221,                0, 0.31981657662220,                0,-0.03407545836221,                0;...
               0,                0,-0.03407545836221,                0,-0.03407545836221,                0,                0;...
               0,                0,                0,-0.00166566566045,                0,                0,                0];


eigd = sqrt(diag(d));
N = size(eigd,1);
eigd = eigd .* randn(N,1);
M = sqrt(N);
phase = reshape(u*eigd,M,M);


rhat = 1/(M-1); % spacing between the samples
%
%
total_points = M;
for i = 1:num_increases
%
% First interpolation
%
   clear filled
   [rows,cols] = size(phase);
   filled (1:2:2*rows-1,1:2:2*cols-1) = phase;
   new = conv2(filled, inter1,'same');

%
%  Next addition of random residual
%
   mask = zeros(2*rows-1,2*cols-1);
   mask(2:2:2*rows-1,2:2:2*cols-1) = ones(rows-1,cols-1);
   new1 = new + mask .* sqrt(residual*rhat^(5/3)) .* ...
               randn(2*rows-1,2*cols-1);

%
% Second interpolation
%
   mask1a = zeros(2*rows-1, 2*cols-1);
   mask1b = zeros(2*rows-1, 2*cols-1);

   mask1a(1:2:2*rows-1, 2:2:2*cols-1) = ones(rows, cols-1);
   mask1b(2:2:2*rows-1,1:2:2*cols-1)  = ones(rows-1,cols);
   mask1 = mask1a + mask1b;
%
   new2 =  conv2(new1, inter2,'same');

   rhat = rhat/sqrt(2);
   phase =  new2 + mask1.*  sqrt(residual*rhat^(5/3)) .* ...
               randn(2*rows-1,2*cols-1);
%
% Calculation to work out the size of the new phase screen after
% truncation to remove edge effects. We llose 5 rows each side.
% New span of phase screen is (2*rows-12)/(2*rows-2)
%
       rhat = rhat/sqrt(2);

       phase = phase(6:2*rows-6,6:2*rows-6);
       size_phase = size(phase,1);

       total_points = 2*total_points - 1;
end

try
scale = (total_points/(size_phase -1))^(5/6);
phase = phase*scale;
catch
    3
end