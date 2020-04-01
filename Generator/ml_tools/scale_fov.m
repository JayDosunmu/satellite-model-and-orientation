function gg = increase_fov( g, scale )

gg = zeros(size(g));
for k=1:size(g,3)
    gg(:,:,k) = p_ImgScale(g(:,:,k), scale);
end

