function compile_evrything_chris(target_folder)



cd(target_folder);

files   = subdir(fullfile(pwd, 'inference.cpp'));
n_files = length(files);

for i = 1 : n_files
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    system('g++ main.cpp optimization.cpp inference.cpp -o InferenceMAP.o -lm');

end
cd(target_folder);


end