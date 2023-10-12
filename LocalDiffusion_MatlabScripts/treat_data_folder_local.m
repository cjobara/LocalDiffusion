function treat_data_folder_local(target_folder, varargin)

%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
prepare_data_REN(target_folder);
cd(target_folder);
%%%%%%%%%%%%%
%%%%%%%%%%%%%
files   = subdir(fullfile(pwd, 'data.mat'));
n_files = length(files);

for i = 1 : n_files
    fprintf('%i\t %i\n', i, n_files);
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    try
    load(files(i).name)
    apply_direct_Mapping_w_tracking_christopher(traj,name_original_file);
    clear traj name_original_file;
    end
    cd(target_folder);
end


%%%%%%%%%%%%
%%%%%%%%%%%%




end