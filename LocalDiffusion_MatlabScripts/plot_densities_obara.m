%function plot_diffusion_map_voronoi(Maps, estimator)
function figure1= plot_densities_obara(Maps, C, varargin)


try
figure1 = figure('XVisual',    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',    'PaperType','a4letter',  'PaperSize',[21 30]);
axes1 = axes('Parent',figure1,'YScale','linear','YMinorTick','off','XScale','linear', 'FontSize',20,'LineWidth',2);
 box(axes1,'on');
hold(axes1,'all');
catch 
    figure1 = figure;
end

if isempty(Maps)
    text(0.25,0.5,'Maps is empty','Color','k','FontSize',30);
else


nombre_couleur = 1e6;
%cmap           = colormap(jet(nombre_couleur+1));

pcolor(Maps.XX, Maps.YY, Maps.densities); colorbar; shading interp;

d_min = min(Maps.densities(:));
d_max = max(Maps.densities(:));

caxis([d_min d_max]);


% for i  = 1 : length(Maps)
% %plot(Maps(i).center_x,Maps(i).center_y,'.k', 'Markersize',15 );
% hold on;
% %fill(Maps(i).voronoi_x, Maps(i).voronoi_y,cmap(valeur(i), :));
% fill(Maps(i).voronoi_x, Maps(i).voronoi_y,D(i) );
% 
% end

if nargin>1
    caxis([C(1) C(2)]);
end


axis([min(Maps.XX(:)) max(Maps.XX(:)) min(Maps.YY(:)) max(Maps.YY(:))]);
h=colorbar;
set(h,'fontsize',20);


ylabel(h, 'Diffusion [\mu m^{2}.s^{-1}]');

% Create xlabel
xlabel('Position [\mu m]','FontSize',24,'FontName','Arial');


% Create ylabel
ylabel('Position [\mu m]','FontSize',24,'FontName','Arial');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












