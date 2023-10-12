function extract_value_ROI_and_around(target_folder, varargin)

warning('off');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0
    target_folder = pwd;
end
cd(target_folder);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files   = subdir(fullfile(pwd, 'point_*'));
n_files = length(files);

for i = 54 : n_files
    fprintf('%i\t %i\n', i, n_files);
    kkk    = strfind(files(i).name, '/'); 
    cd(files(i).name(1:kkk(end)-1) );
    
    %%%%%%%%%%%%%   Load files   %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(files(i).name);
    files2   = subdir(fullfile(pwd, 'Stack_*'));
    load(files2(1).name);
    files2   = subdir(fullfile(pwd, 'Maps_*'));
    load(files2(1).name);
    files2   = subdir(fullfile(pwd, 'data*'));
    load(files2(1).name);

    n_roi = length(polygon_roi);
   
    x_data = tout(:,1);
    y_data = tout(:,2);
    r_data = [x_data, y_data];
    D      = tout(:,6);
    n_D    = length(D);
% 
%     distance_matrix = pdist(r_data);
%     distance_matrix = squareform(distance_matrix);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% in the target roi %%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global_in              = zeros(n_D,1);
    
    for j = 1 : n_roi
        x_poly       = polygon_roi(j).x;
        y_poly       = polygon_roi(j).y;
        in        = inpolygon(x_data,y_data,x_poly,y_poly);
        global_in =  global_in  | in; 
        sum_in    = sum(in);
        if sum(in) == 0
            polygon_roi(j).D_mean   = nan;
            polygon_roi(j).D_median = nan;
            polygon_roi(j).D_std    = nan;
            polygon_roi(j).nb_point = 0;
        else
            polygon_roi(j).D_mean   = mean(D(in));
            polygon_roi(j).D_median = median(D(in));
            polygon_roi(j).D_std    = std(D(in));    
            polygon_roi(j).D_min    = min(D(in));
            polygon_roi(j).D_max    = max(D(in));
            polygon_roi(j).nb_point = sum_in;     
            
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% values in the vincinity %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1 : n_roi
        
        x_poly        = polygon_roi(j).x;
        y_poly        = polygon_roi(j).y;
        x_neigh       = polygon_roi(j).x_neighbours;
        y_neigh       = polygon_roi(j).y_neighbours;
        
        
        
        in_neigh      = inpolygon( x_data , y_data , x_neigh , y_neigh);
        in            = inpolygon( x_data , y_data , x_poly  , y_poly );
        
        in_neigh_pure = (in_neigh & ~in );
        sum_in_neigh_pure       = sum(in_neigh_pure);
        
        polygon_roi(j).D_neigh_mean   = mean(D(in_neigh_pure));
        polygon_roi(j).D_neigh_median = median(D(in_neigh_pure));
        polygon_roi(j).D_neigh_std    = std(D(in_neigh_pure));    
        polygon_roi(j).D_neigh_min    = min(D(in_neigh_pure));
        polygon_roi(j).D_neigh_max    = max(D(in_neigh_pure));
        polygon_roi(j).nb_neigh_point = sum_in_neigh_pure;   
        
        
        plot_associated_results(Maps, X, Y, Image_mean_normalized, polygon_roi,j  );
        
        
    end    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% in the intermediate target roi %%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    global_intermediate_in = zeros(n_D,1);
    n_inter_roi            = length(polygon_intermediate_roi);
    
    for j = 1 : n_inter_roi
        x_poly                  = polygon_intermediate_roi(j).x;
        y_poly                  = polygon_intermediate_roi(j).y;
        in                      = inpolygon(x_data,y_data,x_poly,y_poly);
        global_intermediate_in  =  global_intermediate_in  | in; 
        sum_in                  = sum(in);
        if sum(in) == 0
            polygon_intermediate_roi(j).D_mean   = nan;
            polygon_intermediate_roi(j).D_median = nan;
            polygon_intermediate_roi(j).D_std    = nan;
            polygon_intermediate_roi(j).D_min    = nan;
            polygon_intermediate_roi(j).D_max    = nan;       
            polygon_intermediate_roi(j).nb_point = 0;
            
        else
            polygon_intermediate_roi(j).D_mean   = mean(D(in));
            polygon_intermediate_roi(j).D_median = median(D(in));
            polygon_intermediate_roi(j).D_std    = std(D(in));     
            polygon_intermediate_roi(j).D_min    = min(D(in));
            polygon_intermediate_roi(j).D_max    = max(D(in)); 
            polygon_intermediate_roi(j).nb_point = sum_in;       
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% outside the roi %%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    global_not_in   = ~( global_in | global_intermediate_in  );
    sum_not_in      = sum(global_not_in);
    D_not_in_mean   = mean(D(global_not_in));
    D_not_in_median = median(D(global_not_in));
    D_not_in_std    = std(D(global_not_in));
    D_not_in_min    = min(D(global_not_in));
    D_not_in_max    = max(D(global_not_in));
    
    for j = 1 : n_roi
        polygon_roi(j).D_not_in_mean    = D_not_in_mean;
        polygon_roi(j).D_not_in_median  = D_not_in_median;
        polygon_roi(j).D_not_in_std     = D_not_in_std;
        polygon_roi(j).D_not_in_min     = D_not_in_min;
        polygon_roi(j).D_not_in_max     = D_not_in_max;
        polygon_roi(j).nb_point_outside = sum_not_in ;  
        
    end
    
 
%     for i 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save(['results_' files(i).name(kkk(end)+1: end-4) '.mat'], 'polygon_roi', 'polygon_intermediate_roi','Maps', 'X', 'Y', 'Image_mean_normalized','x_data','y_data');
    clear polygon_roi polygon_intermediate_roi Maps X Y Image_mean_normalized x_data y_data;
    cd(target_folder);
    
    
end


end