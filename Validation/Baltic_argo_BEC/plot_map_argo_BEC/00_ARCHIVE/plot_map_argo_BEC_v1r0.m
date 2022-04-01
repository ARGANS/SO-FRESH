% Syntax: plot_map_argo_BEC (script-snipt of Baltic_argo_BEC_analyses.m)
%
% Description
% Plot location of BEC-ARGO floats within a given region and zoom-in to the
% flaot location.
% 
% Version: v1r0
% ========================================================================



% Plots in one map locations of same platform

%% [1] Plot Argo float location within region
if nn > 1 && platform(nn) ~= platform(nn-1)
    figure(2); hold on
else
    figure(2); clf; hold on
end

set(gcf,'DefaultAxesFontSize',24);
% set(gca,'TickLabelInterpreter','tex');


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


title({['Argo ' basin_str ' ('...
    num2str(lat_min) '-' ...
    num2str(lat_max) '\circ N) ' ]; ...
    [datestr(itime_start,'yyyymmdd') '-'...
    datestr(itime_end,'yyyymmdd')]});
hold on

h1 = plotm(lat,lon,'ko','MarkerFaceColor','r',...
    'MarkerSize',12); hold on

h2 = plotm(lat(nn),lon(nn),...
    'ks','linewidth',1,'MarkerSize',10,...
    'MarkerFaceColor','b');


% Save figure  - output -
fg_name = [fn '_MAP_PLAT'...
    num2str(platform(nn)) ...
    '_PROF#' num2str(ID(nn)) '.'  fg_format];


% Save figure  - output -
folder_this = [folder_figs 'MAPS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];


% check fn existence
fg_exist = exist(fg_name,'file');
if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end

%% [2] Plot ZOOM-IN to Argo float location
if nn > 1 && platform(nn) ~= platform(nn-1)
    figure(22); hold on
else
    figure(22); clf; hold on
end

set(gcf,'DefaultAxesFontSize',24);


title({['Argo ' basin_str ' ('...
    num2str(lat_min) '-' ...
    num2str(lat_max) '\circ N) ' ]; ...
    [datestr(itime_start,'yyyymmdd') '-'...
    datestr(itime_end,'yyyymmdd')]});
hold on

h1 = plot(lat,lon,'ko','MarkerFaceColor','r',...
    'MarkerSize',12); hold on

h2 = plot(lat(nn),lon(nn),...
    'ks','linewidth',1,'MarkerSize',10,...
    'MarkerFaceColor','b');
fillmap

grid on
box on


% Save figure  - output -
fg_name = [fn '_MAP-ZOOM_PLAT'...
    num2str(platform(nn)) ...
    '_PROF#' num2str(ID(nn)) '.'  fg_format];


% Save figure  - output -
folder_this = [folder_figs 'MAPS/MAP-ZOOM/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];


% check fn existence
fg_exist = exist(fg_name,'file');
if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end





