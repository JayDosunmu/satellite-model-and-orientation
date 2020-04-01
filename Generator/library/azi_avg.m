function [ r_avg] = azi_avg( frame,Fmask ) 

% function [ r_avg] = azi_avg( xxframe ) 
% return the radial averages for an image or something square...

xdim = size(frame,1);   %** assume square image
start = xdim/2;

[x,y]   = meshgrid(-start:1:start-1,-start:1:start-1);
r_dist  = sqrt(x.^2 + y.^2);
rad_avg = zeros(start,1);
t       = r_dist(find(r_dist > 0)); 
max_rad = round(max(max(r_dist)))+1;
t_bin   = hist(t,max_rad);
r_avg   = zeros(start,1);
for i=0:start-1
    idx = find(round(r_dist) == i & Fmask ==1);
    r_avg(i+1) = sum(frame(idx))./length(idx);
end
    
