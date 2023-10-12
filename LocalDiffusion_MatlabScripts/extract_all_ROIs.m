function extract_all_ROIs(target_folder,criterion_roi, varargin)
warning('off')

%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
%%%%%%%%%%%%%
%%%%%%%%%%%%%
files   = subdir(fullfile(pwd, 'Stack_*'));
n_files = length(files);


 for i = 1 : n_files 
    try
        load(files(i).name);
        
        
        
        
    catch
       fprintf('oops\n') 
    end
 end



end