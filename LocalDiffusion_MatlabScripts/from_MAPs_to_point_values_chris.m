function from_MAPs_to_point_values_chris(target_folder, varargin)



cd(target_folder);

files   = subdir(fullfile(pwd, 'Maps*'));
n_files = length(files);

for i = 1 : n_files
    
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    name_root = files(i).name(kkk(end)+1:end-4);
    
    try
        load(files(i).name);
        n_Maps = length(Maps);
        fichier = fopen(['position_value_' name_root '.txt'], 'w+');
        for j = 1 : n_Maps;
            n_x = length(Maps(j).x);
            for k = 1 : n_x
                fprintf(fichier, '%f\t %f\t %f\t %f\t', Maps(j).x(k), Maps(j).y(k), Maps(j).t(k) ,Maps(j).D);
                try
                    fprintf(fichier, '%f\t %f\n', Maps(j).Fx, Maps(j).Fy);
                catch
                    fprintf(fichier,'\n');
                end
                
            end
        end
        
        
    end
end

cd(target_folder);


end