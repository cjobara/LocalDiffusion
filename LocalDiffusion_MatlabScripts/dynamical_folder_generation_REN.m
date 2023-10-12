function [name_folder, filename,pathname, inum] =  dynamical_folder_generation_REN(i, files)

inum = num2str(i);
name   = files(i).name;


if ismac||isunix
    kkk    = strfind(name, '/'); 
elseif ispc
    kkk    = strfind(name, '\');  
end

cd(files(i).name(1:kkk(end)-1) );

name_folder = ['folder_' name(kkk(end)+1:end-4)];
mkdir(name_folder);

filename = name(kkk(end):end);
pathname = name(1:kkk(end)-1);

% fprintf('path name %s\n',pathname );
% fprintf('file name %s\n',filename );

end