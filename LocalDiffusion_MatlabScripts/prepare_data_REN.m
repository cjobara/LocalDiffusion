function prepare_data_REN(target_folder, varargin)

%%
if nargin ==0
target_folder     = pwd;
end
cd(target_folder);
%%
files  = subdir(fullfile(pwd, '*.txt'));
n_files = length(files);
%%
for i =  1 : n_files
    [name_folder, filename,pathname, inum] =  dynamical_folder_generation_REN(i, files);
    cd(name_folder);

%     [dt_theo,nb, t, x, y]                  = import_from_sengupta_files(files(i).name);
%     [dt_theo,nb, t, x, y]                  = import_from_christopher_files(files(i).name);
try
    [dt_theo,nb, t, x, y]                  = import_from_christopher_files_v2(files(i).name);
    t                                      = t*dt_theo(1);
    traj                                   = [nb, x, y, t];
    name_original_file                     = filename(2:end);
    save('data.mat', 'traj', 'name_original_file');
end
    cd(target_folder);
end


cd(target_folder);




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












