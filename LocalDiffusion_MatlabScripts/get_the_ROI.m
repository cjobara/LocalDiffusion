function [polygon_roi,polygon_intermediate_roi] = get_the_ROI(Image_mean_normalized,X, Y, criterion_roi)

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

critere_roi_low = 0.25;
nb_pixel        = 128;
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%


sigma_pixel = 3;
densities   = Image_mean_normalized;
G           = fspecial('gaussian',[10 10],sigma_pixel);
densities   = imfilter(densities,G,'same');


density_criterion = criterion_roi;


[cent]              = FastPeakFind(densities, density_criterion );
intensities = [];
for j = 1 : floor(length(cent)/2)
    intensities = [intensities; densities(cent(2*j),cent(2*j-1))];
end




indice = 1;
for j = 1 : floor(length(cent)/2)
    if (intensities(j)>= density_criterion)
    x_target(indice,1) =   X(cent(2*j), cent(2*j-1));
    y_target(indice,1) =   Y(cent(2*j), cent(2*j-1));
    indice             = indice + 1;
    end
end

Image_mean_normalized_th      = Image_mean_normalized;
II                            = Image_mean_normalized> criterion_roi;
Image_mean_normalized_th(II)  = 1; 
Image_mean_normalized_th(~II) = 0;
CC                            = bwconncomp(Image_mean_normalized_th);
n_roi                         = CC.NumObjects;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B                             = boundaries(Image_mean_normalized_th, 8);
n_roi                         = length(B);
n_roi_domains                 = n_roi;

for i = 1 : n_roi 
    x_roi                   = (B{i,1}(:,2)-1)*0.16;
    y_roi                   = (B{i,1}(:,1)-1)*0.16;
    polygon_roi(i).x        = x_roi;
    polygon_roi(i).y        = y_roi;
    A                       = polyarea(x_roi,y_roi);
    polygon_roi(i).surface  = A;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Image_mean_normalized_th      = Image_mean_normalized;
II                            = (Image_mean_normalized<= criterion_roi) & (Image_mean_normalized>critere_roi_low);
Image_mean_normalized_th(II)  = 1; 
Image_mean_normalized_th(~II) = 0;

B                             = boundaries(Image_mean_normalized_th, 8);
n_roi                         = length(B);

for i = 1 : n_roi 
    x_roi                                = (B{i,1}(:,2)-1)*0.16;
    y_roi                                = (B{i,1}(:,1)-1)*0.16;
    polygon_intermediate_roi(i).x        = x_roi;
    polygon_intermediate_roi(i).y        = y_roi;
    A                                    = polyarea(x_roi,y_roi);
    polygon_intermediate_roi(i).surface  = A;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1 : n_roi_domains 
    
    BW                            = zeros(nb_pixel,nb_pixel);
    pixel_roi_x                   = floor(polygon_roi(i).y/0.16) + 1;
    pixel_roi_y                   = floor(polygon_roi(i).x/0.16) + 1;
    for j = 1 : length(pixel_roi_x);
        BW(pixel_roi_x(j), pixel_roi_y(j)) = 1;
    end
    clear pixel_roi_x pixel_roi_y
    
    BW2 = imfill(BW,'holes');
    SE2 = strel('diamond',10);
    BW3 = imdilate(BW2,SE2);
    B                             = boundaries(BW3, 8);
    
    x_roi_neighbours               = (B{1,1}(:,2)-1)*0.16;
    y_roi_neighbours               = (B{1,1}(:,1)-1)*0.16;
    polygon_roi(i).x_neighbours    = x_roi_neighbours;
    polygon_roi(i).y_neighbours    = y_roi_neighbours;
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%     K                = convhull(xx, yy);
%     x_roi            = xx(K);
%     y_roi            = yy(K);
%     polygon_roi(i).x = x_roi;
%     polygon_roi(i).y = y_roi;
