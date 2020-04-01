function [g4, v_az, v_el] = gen_object_poses(stl_fname, n_poses, dsf, v_az, v_el ) 

% g = gen_object_poses(stl_fname, n_poses, dsf,  v_az, v_el ) 
%
% stl_fname         file name to STL file (full path)
% n_poses           number of random poses
% dsf               downsampling factor
% v_az              vector azimuthal camera positions
% v_el              vector of elevation camera positions
%
% returns: g        NxNxn_poses 3D array of object poses
%

if nargin==3
    % Camera view angle (degrees) random
    v_az =  rand(1,n_poses).*180;
    v_el =  rand(1,n_poses).*180;
end
n_poses = length(v_az);


% Camera view angle (degrees) sequence
% v_az = -30:1:0;
% v_el = 40:1:70;

%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure

% stl_fname = '/Users/douglashope/Dropbox/hsr/ML/Objects STL/cassini.stl';
fits_fname = strrep(stl_fname,'.stl','.fits');

fv = stlread1(stl_fname);


%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.
figure(5)
clf
set(gcf,'Color','w');
patch('Faces',fv.faces,'Vertices',fv.vertices,'FaceColor',[0.8 0.8 1.0], 'EdgeColor','none','FaceLighting', 'gouraud','AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('left');
material('dull');
lighting('GOURAUD');
% Fix the axes scaling, and set a nice view angle
axis('image');axis off

set(gca,'CameraTarget',[0 0 0]);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
camtarget('manual')

view([v_az(1) v_el(1)]);
gg=zeros(896,896,length(v_az));

parfor k=1:length(v_az)
    az_view = v_az(k);
    el_view = v_el(k);
    view(az_view,el_view);
    camva(10)
        
    FF=getframe;
    g1 = double(rgb2gray(FF.cdata));
    g1 = abs(g1 - g1(1,1));

    pad_size = (896 - size(g1))/2;
    
    if pad_size(1) ~= fix(pad_size(1))
        pad_size1_pre  = fix(pad_size(1));
        pad_size1_post    = 896-size(g1,1)-pad_size1_pre;
    else
         pad_size1_pre  = pad_size(1);
         pad_size1_post = pad_size(1);
    end
    if pad_size(2) ~= fix(pad_size(2))
        pad_size2_pre  = fix(pad_size(2));
        pad_size2_post    = 896-size(g1,2)-pad_size2_pre;
    else
        pad_size2_pre  = pad_size(2);
        pad_size2_post = pad_size(2);
    end
    pad_pre     = [ pad_size1_pre pad_size2_pre ];
    pad_post    = [ pad_size1_post pad_size2_post ];
    g2          = padarray(padarray(g1, pad_pre, 'pre'), pad_post, 'post');
    gg(:,:,k) = g2;
end

if dsf>1 
    g4=sepblockfun(gg,[dsf dsf],@mean);
    h=fspecial('Gaussian',1);
    g4 = imfilter(g4,h);
else
    g4 = gg;
end