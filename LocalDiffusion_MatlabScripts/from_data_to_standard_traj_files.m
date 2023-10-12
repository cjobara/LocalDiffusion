function from_data_to_standard_traj_files()


cd(target_folder);

files   = subdir(fullfile(pwd, 'data.mat'));
n_files = length(files);

for i = 1 : n_files
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    try
    load(files(i).name)
    name_file = files(i).name(kkk(end)+1:end);
    plot_diffusion_maps_voronoi_onto_points_to_figure(Maps);
    saveas(gcf,[name_file(1:end-4) '.eps'], 'epsc') ;
    close all force;
    end
end






end