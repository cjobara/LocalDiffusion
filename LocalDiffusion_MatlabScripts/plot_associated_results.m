function plot_associated_results(Maps, X, Y, Image_mean_normalized, polygon_roi,j  )

% try
% figure1 = figure('XVisual',    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',    'PaperType','a4letter',  'PaperSize',[21 30]);
% axes1 = axes('Parent',figure1,'YScale','linear','YMinorTick','off','XScale','linear', 'FontSize',20,'LineWidth',2);
%  box(axes1,'on');
% hold(axes1,'all');
% catch 
%     figure1 = figure;
% end

j_num = num2str(j);

figure1 = plot_diffusion_maps_voronoi_onto_points(Maps);
hold on;
% pcolor(X,Y, Image_mean_normalized); shading interp; colorbar;
hold on;
plot(polygon_roi(j).x,polygon_roi(j).y,'.-k', 'LineWidth', 4 );
plot(polygon_roi(j).x_neighbours,polygon_roi(j).y_neighbours,'.-m', 'LineWidth', 4 );
hold on ;



h=colorbar;
set(h,'fontsize',20);

xlabel({'Position [\mu m]'})
ylabel({'Position [\mu m]'})


min_x = min(min(X));
max_x = max(max(X));

min_y = min(min(Y));
max_y = max(max(Y));


axis([min_x max_x min_y max_y]);

saveas(figure1,['diffusion_values_domain_roi_' j_num '.eps'], 'epsc') ;
close all force;



end