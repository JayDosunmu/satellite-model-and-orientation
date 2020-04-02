addpath('library');
addpath('ml_tools');
addpath('asim');
addpath('objects_training')

% Where are the STL files stored 
object_dir = 'objects_training/';

%object_dir = fullfile('E:\','objects_training',filesep);
Dr0        = 5;            % turbulence strength D/r0 =5 (low) /D/r0 = 15 (medium) / D/r0=21 (strong)
n_poses    = 100;             % how many random or specified poses 
dsf        = 4;             % downsample factor (controls image grid size)
fov        = 0.5;           % field of view

% generate specific object poses 
% v_az = zeros(1,11);
% v_el = 0:10:100;
% n_poses = length(v_el);

% measure elapsed time
t = cputime;

adir = dir([ object_dir '*.stl']);
for ss=1:length(adir)
    if adir(ss).isdir==0
        [toss,name,ext] = fileparts(adir(ss).name);
        name = strrep(name, "._", "");
        stl_fname   = join([ object_dir filesep name ext ], '');
        
       [g, v_az, v_el]   = gen_object_poses(stl_fname, n_poses, dsf);  % generate random poses 
%        [g, ~, ~]  = gen_object_poses(stl_fname, n_poses, dsf, v_az, v_el);  % generate specific poses 
        g2   = scale_fov(g, 0.5);           % change the FOV - make object smaller<1 /bigger >1
        gp  = blur_object_poses(g2, Dr0);       % Specify turbulence strength via D/r0 value
        
        for n=1:n_poses
            imwrite(gp(:,:,n), join(["output" filesep name "_" int2str(n) ".JPEG"], ''), "JPEG");
        end
    end
end

% record elapsed time
e = cputime -t;

% diplay elapsed time
disp('elapsed time:')
disp(e)
