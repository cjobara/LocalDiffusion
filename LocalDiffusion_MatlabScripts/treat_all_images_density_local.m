function treat_all_images_density_local(target_folder, varargin)

%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
%%%%%%%%%%%%%
%%%%%%%%%%%%%
files   = subdir(fullfile(pwd, 'densities_*'));
n_files = length(files);

for i = 1 : n_files
    
    fprintf('%i\t %i\n ', i, n_files)
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    try
        load(files(i).name)
        files2   = subdir(fullfile(pwd, 'densities_*'));
        name_loc = files2(1).name(1:end-4);
    %     apply_direct_Mapping_w_tracking_christopher(traj,name_original_file);
        h = plot_densities_obara(Maps_densities);  
        saveas(h,[name_loc '.eps'], 'epsc') ;
        close all force;
        clear Maps;
    end
    cd(target_folder);
    
end


%%%%%%%%%%%%
%%%%%%%%%%%%




end