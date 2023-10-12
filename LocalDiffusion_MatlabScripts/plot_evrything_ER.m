function plot_evrything_ER(target_folder)



cd(target_folder);

files   = subdir(fullfile(pwd, 'Maps_*'));
n_files = length(files);

for i = 1 : n_files
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    try
    load(files(i).name)
    name_file = files(i).name(kkk(end)+1:end);
    plot_diffusion_maps_voronoi_onto_points_to_figure(Maps);
    saveas(gcf,['diffusion_' name_file(1:end-4) '.eps'], 'epsc') ;
    close all force;
    end
end


cd(target_folder);

files   = subdir(fullfile(pwd, 'Maps_*'));
n_files = length(files);

for i = 1 : n_files
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    try
    load(files(i).name)
    name_file = files(i).name(kkk(end)+1:end);
    plot_norm_force_maps_voronoi_onto_points_to_figure(Maps);
    saveas(gcf,['velocity_' name_file(1:end-4) '.eps'], 'epsc') ;
    close all force;
    end
end




end