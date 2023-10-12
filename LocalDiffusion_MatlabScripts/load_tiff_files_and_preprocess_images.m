function load_tiff_files_and_preprocess_images(target_folder, varargin)
warning('off')



criterion_roi = 0.9;

%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
%%%%%%%%%%%%%
%%%%%%%%%%%%%
files   = subdir(fullfile(pwd, '*.tif'));
n_files = length(files);



for i = 1 : n_files
    try
    fprintf('%i\t %i\n', i, n_files);
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
%     load(files(i).name)
    
    FileTif               = files(i).name;
    InfoImage             = imfinfo(FileTif);
    mImage                = InfoImage(1).Width ;
    nImage                = InfoImage(1).Height ;
    NumberImages          = length(InfoImage) ;
    FinalImage            = zeros(nImage,mImage,NumberImages,'uint8');
    for j = 1 : NumberImages
        FinalImage(:,:,j) = imread(FileTif,'Index',j);
    end
    Image_mean            = mean(FinalImage, 3);
    Image_mean_normalized = mean(FinalImage, 3)./255;
    xx    =  [0:1:127]*0.16;
    yy    =  [0:1:127]*0.16;
    [X,Y] = meshgrid(xx,yy);
    
    
    [polygon_roi,polygon_intermediate_roi] = get_the_ROI(Image_mean_normalized,X, Y, criterion_roi);
    
    
    
    save(['Stack_' files(i).name(kkk(end)+1: end-4) '.mat'], 'FinalImage','NumberImages', 'nImage', 'mImage', 'InfoImage','FileTif' , ...
        'Image_mean', 'Image_mean_normalized', 'xx', 'yy', 'X', 'Y', 'polygon_roi', 'polygon_intermediate_roi');
    
    load('data.mat');
    
    figure1 = figure('XVisual',    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',    'PaperType','a4letter',  'PaperSize',[21 30]);

    pcolor(X,Y, Image_mean_normalized); shading interp; colorbar;
    hold on ;
    plot(traj(:,2), traj(:,3), '.-k', 'Markersize', 2, 'LineWidth', 1);
    saveas(figure1,['overlay_' name_original_file '.eps'], 'epsc') ;
    close all force;
    close all force;
    pause(1);
    
    clear test FinalImage NumberImages nImage mImage FileTif InfoImage Image_mean Image_mean_normalized ;
    clear traj name_original_file
    
    
    end
    cd(target_folder);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








