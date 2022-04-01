function h1 = plot_map_SDN_BEC(lon,lat)

% Syntax: plot_map_seadatanet_BEC (script-snipt of Baltic_seadatatanet_BEC_analyses.m)
%
% Description
% Plot location of SeaDataNet stations within a given region and zoom-in to the
% flaot location.
% 
% Version: v1r1
% ========================================================================



% Plots in one map locations of same platform

%% [1] Plot Argo float location within region
set(gcf,'DefaultAxesFontSize',24);

[xmin,xmax,ymin,ymax,basin_str] = map_lim_raf (9);



Y1 = lon;
X1 = lat;
% Z1 = time_range(ind_platform);

% itime_start = time_range(1);
% itime_end = time_range(end);
% itime_start_str = datestr(itime_start,'yyyymmdd');
% itime_end_str = datestr(itime_end,'yyyymmdd');


map_projection = 'merc';

lat_min = ymin;
lat_max = ymax;
lon_min = xmin;
lon_max = xmax;

lon_step = 10;
lat_step = 5;

fillmap_super(map_projection,...
    lon_min,lon_max,lat_min,...
    lat_max,lon_step,lat_step);
hold on

h1 = plotm(X1,Y1,...
    'ko','MarkerFaceColor','r','MarkerSize',5); hold on





% %% [2] Plot ZOOM-IN to Argo float location
% if nn > 1 && platform_irange(nn) ~= platform_irange(nn-1)
%     figure(22); hold on
% else
%     figure(22); clf; hold on
% end
% 
% set(gcf,'DefaultAxesFontSize',24);
% 
% 
% title({['Argo ' basin_str ' (ZOOM-IN '...
%     num2str(lat_min) '-' ...
%     num2str(lat_max) '\circ N) ' ]; ...
%     [datestr(itime_start,'yyyymmdd') '-'...
%     datestr(itime_end,'yyyymmdd')]});
% hold on
% 
% h1 = plot(X1,Y1,'ko','MarkerFaceColor','r',...
%     'MarkerSize',12); hold on
% 
% h2 = plot(X1(nn),Y1(nn),...
%     'ks','linewidth',1,'MarkerSize',10,...
%     'MarkerFaceColor','b');
% fillmap
% 
% grid on
% box on
% 
% 
% % Save figure  - output -
% fg_name = ['argo_' datestr(Z1(nn),'yyyymmdd') '_MAP-ZOOM_PLAT'...
%     num2str(platform_irange(nn)) ...
%     '_PROF#' num2str(ID_irange(nn)) '.'  fg_format];
% 
% 
% % Save figure  - output -
% folder_this = [folder_figs 'MAPS/MAP-ZOOM/'];
% 
% if fg_save == 1
%     foldercheck_raf(folder_this); %! make folder_figs
% end
% 
% fg_name = [folder_this fg_name];
% 
% 
% % check fn existence
% fg_exist = exist(fg_name,'file');
% if fg_save == 1 && fg_exist == 0
%     save_raf(gcf,fg_name,fg_format); close
% end
% 
% 



