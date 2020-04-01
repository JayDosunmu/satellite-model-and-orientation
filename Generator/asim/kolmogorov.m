
function [screen] = kolmogorov(grid_size, D_ro);
%  Simulation of a phase screen with Kolmogorov Statistics
%
%  function [screen] = kolmogorov(grid_size, D_ro);
%
%  Inputs       grid_size       - desired phase screen size
%               D_ro            - D/r0 value, determines severity of turbulence
%
%  Outputs      screen          - a phase screen with Kolmogorov Statistics
%                                 of size [grid_size, grid_size]
%
%  Requires     genphase.m
%               svectors4.mat
%
%
%  Rachel Johnston, Cressida Harding and Richard Lane April 1999
%  Department of Electrical and Electronic Engineering
%  University of Canterbury
%  Christchurch
%  New Zealand
%
%  Code may be copied and used freely provided acknowledgement of the original source
%  is maintained.
%
% Finding the number of increases required
% Starting from a phasescreen of 15 by 15
start_size  = 15;
no_increases = 0;
final_size = 0;
while (final_size < grid_size)
       final_size = 2*start_size - 11;
       start_size = final_size;
       no_increases = no_increases + 1;
end

if no_increases == 0 
   3 
end
%Generating the phasescreen
full_screen = genphase(no_increases);
chopped_screen = full_screen(1:grid_size, 1:grid_size);

%Need to rescale the screen so that it has the appropriate D_ro after
%the chopping
size_chopped = size(chopped_screen,1);
size_full = size(full_screen,1);
scale = (size_full/(size_chopped - 1))^(5/6);
scaled_screen = chopped_screen*scale;

%Scaling so it has the specfied D_ro
for b=1:length(D_ro)
    screen(:,:,b) = scaled_screen*((D_ro(b))^(5/6));
end
