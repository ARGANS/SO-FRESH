% PLOT: Grid Argo locations to the SMOS-GRID
%
%
% Syntax: Baltic_argo_BEC_colocations_PLOT_griddedARGO.m
%
% =========================================================================


figure

% [1] Original SDN
% h1  = plot(lon(nn),lat(nn),'ko','MarkerFaceColor','r','MarkerSize',5); 
h1  = plot(lon_alfa,lat_alfa,'ko','MarkerFaceColor','r','MarkerSize',5); 

hold on


% use r2 as the modified searching radius (+25% of r) at high-lats
% circ_01 = plot_study_region(lon(nn),lat(nn),r2);
circ_01 = plot_study_region(lon_alfa,lat_alfa,r2);
set(circ_01,'color','r','linewidth',1.5)

% [2] Argo float on smos-grid
% h2 = plot(lon_grid(nn),lat_grid(nn),'sb','Markersize',15,'linewidth',1);
h2 = plot(lon_grid,lat_grid,'sb','Markersize',15,'linewidth',1);


% circ_02 = plot_study_region(lon_grid(nn),lat_grid(nn),r2);
circ_02 = plot_study_region(lon_grid,lat_grid,r2);

set(circ_02,'color','b','linewidth',1.5)

% [3] smos-grid
h3 = plot(lon_sss,lat_sss,'+','color',[1 1 1].*.5);
hold on

xlim([5 xmax]); ylim([ymin 68]);

fillmap

lg = legend([h1(1),h2(1),h3(1)],{'orginal','gridded','smos-grid'},...
    'location','Southeast','fontsize',18);



title('Baltic REF-2-SAT grid','fontsize',24)
xlabel ('longitude (\circ)','fontsize',24);
ylabel ('latitude (\circ)','fontsize',24);