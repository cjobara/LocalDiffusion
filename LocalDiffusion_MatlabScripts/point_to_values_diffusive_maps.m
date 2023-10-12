function point_to_values_diffusive_maps(target_folder, varargin)
warning('off');


%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
%%%%%%%%%%%%%
%%%%%%%%%%%%%
files   = subdir(fullfile(pwd, 'Maps_*'));
n_files = length(files);

for i = 1 : n_files
    fprintf('%i\t %i\n', i, n_files);
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    
    load(files(i).name);
    tout = [];
    for j = 1 : length(Maps)
        x           = Maps(j).x;
        n_x         = length(x);
        y           = Maps(j).y;
        t           = Maps(j).t;
        center_x    = Maps(j).center_x*ones(n_x,1);
        center_y    = Maps(j).center_y*ones(n_x,1);
        D_vec       = Maps(j).D*ones(n_x,1);
        tout        = [tout; x, y, t, center_x,  center_y,D_vec ];
        
    end
    
    
    save(['point_values_' files(i).name(kkk(end)+1: end-4) '.mat'], 'tout');
    clear x n_x y t center_x center_y D_vec tout;
        
end



cd(target_folder);

end
