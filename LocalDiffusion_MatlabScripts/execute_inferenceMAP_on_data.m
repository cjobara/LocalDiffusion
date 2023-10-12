function execute_inferenceMAP_on_data(target_folder, varargin)



cd(target_folder);

files   = subdir(fullfile(pwd, '*.cluster'));
n_files = length(files);


for i  = 1 : n_files
    
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    name_file = files(i).name(kkk(end)+1:end-8);
    system(['./InferenceMAP.o -i ' name_file '.cluster -o ' name_file '.vmesh -e 20 -m Dvel -d 0. -v 0. -j 0 ']);
    
end
cd(target_folder);


end