function apply_direct_Mapping_w_tracking_christopher(data,name_original_file)

%% local parameters
d        =  2;
D_thresh =  3.;


traj = data;
name_input_file = name_original_file(1:end-4);
full_name_input1 = [name_input_file '.txt'];

fprintf('traj loaded\n');

%% parameters general 
 name_pipeline = 'global_maps';
 initial_mesh = [];
 Previous_Assignment = [];
 state             = 'diff_voronoi';
 modal             = 'hungarian';
 criterion         = 'number_point_sliding';
%% parameters and update to do in the generate_parameters 
%%
t_loc     = traj(2:end,4) - traj(1:end-1,4);
II        = t_loc <= 0;
t_loc(II) = nan;
dt_theo   = nanmin(t_loc);

% fprintf('%f\n',dt_frame);
%% take_fulll_trajectory
%%
[xyt,t_min, t_max] = generate_xyt_neuron(traj);
movie_per_frame    = convert_traj_struct(xyt);
movie_per_frame    = regularize_movie_per_frame_clusters_VLP(movie_per_frame, t_min, t_max, dt_theo);

fprintf('movie per frame\n');
%%
%%
parameters        = generate_parameters(movie_per_frame, [], name_pipeline);
%%
%%

number_per_zone_num     = num2str(parameters.number_per_zone);
number_min_per_zone_num = num2str(parameters.min_number_per_zone);
name_output_file        = [name_input_file '_kmeans_' number_per_zone_num '_nb_min_' number_min_per_zone_num];
name_Maps_output        = ['Maps_' name_output_file '.mat'];
name_densities_output   = ['densities_' name_output_file '.mat'];
fprintf('output\n');
%%
%%
tout                   = generate_tout_from_predefine_trajectories(traj, d);
tout                   = clean_tout_threshold(tout, D_thresh, dt_theo,d);

fprintf('tout\n');
[Maps, Maps_densities] = Mapping_without_tracking_from_tout(tout,state,parameters, dt_theo);
fprintf('maps\n');
generate_file_for_pure_optimization(Maps, name_output_file);
save(name_Maps_output, 'Maps');
save(name_densities_output, 'Maps_densities');

end


