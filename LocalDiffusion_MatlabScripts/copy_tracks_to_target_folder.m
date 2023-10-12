function copy_tracks_to_target_folder(target_folder, varargin)

%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
%%%%%%%%%%%%%
%%%%%%%%%%%%%


files   = subdir(fullfile(pwd, '*Tracks.txt'));
n_files = length(files);


for i=1 : n_files
   name_loc = files(i).name;
   kkk    = strfind(files(i).name, '/'); 
   cd(files(i).name(1:kkk(end)-1) );
   
   name_without_extension = name_loc( kkk(end)+1  :end-4);
   name_with_extension     = name_loc( kkk(end)+1  :end);
   name_tracks            = 'Tracks';
   
   
   fprintf('%s\n', name_without_extension);
   fprintf('%s\n', ['mv ' name_with_extension  '  folder_' name_without_extension(1:end-6) name_tracks '/' name_with_extension  ]);
   
   system(['mv ' name_with_extension  '  folder_' name_without_extension(1:end-6) name_tracks '/' name_with_extension  ]);
%     fprintf('%s\n',['mv ' name_with_etension  '  folder_' name_without_extension(1:end-8) name_tracks '/' name_with_etension  ] );
end
    
end