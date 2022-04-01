function Baltic_plot_fun(ibasin,itime,data_type,varargin)
%
% Syntax: Baltic_plot_fun(ibasin,itime,data_type,plot_type,region)
% 
% Description
% Plot to see SSS BEC products in the Arctic.
% This function builds up on the script Arctic_plot_v1r1.m
%
% Input (as seen in function rd_smos_L4_BEC_v1r2.m)
% itime: date as in MATLAB number
% ibasin: basin number (for more info type help map_lim_raf)
% data_type for Arctic: 
%           [1] 9_days Arctic+ product (v3.0) 
%           [2] 9_days Arctic BEC product (v2.0)
%           [3] Global BEC product (v001)
%           [4] Global (default) 
%
% plot_type for Arctic
%           [1] Mercartor
%           [2] Polar 
%           [3] Lambert [plot each side of the Arctic: Atlantic, Pacific
%           sides]
%           [4] All- Plot all projection above
%
% Version: v1r0
% Author: rcatany (rafaelcatany@gmail.com)
% 
% History
% Version | Date| Note
% v1r0 [24/10/2019] : Creation of the function
% =======================================================================
%%

% clc;
% clear

% Example (settings)
% ibasin = 7;     % Basin number (Arctic, ibasin 7)
% data_type = 3;  % [1-3] Arctic+ ([1]monthly,[2] 7day or [3] 9 day) product and [4] Global product
% 
% iyear = 2014;
% imonth = 2;
% iday = 1;
% 
% itime = datenum(iyear,imonth,iday);

if isempty(varargin)
    plot_type = 2; % Mercator projection *set default
    region = 'ALL';
end

% set the plot type
if ~isempty(varargin{1})
plot_type = varargin{1};
end

% Especify section of the Arctic to focus: [1] NA; [PA]; or [ALL] sections
if ~isempty (varargin{2})
    region = varargin{2};
end


[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

if data_type == 1
    % ARC_PRODUCT = ibasin == 7 &&  data_type == 1;% Arctic+ product (v3.0)
    data_type_str = ['ARCTIC+ product (v3.0)'];
    
elseif data_type == 2
    % ARC_PRODUCT = ibasin == 7 && data_type == 2; % Arctic BEC prodcut (v2.0)
    data_type_str = ['ARCTIC product (BEC v2.0)'];
else
    % ARC_PRODUCT = 0;
    data_type_str = ['GLOBAL product (v001)'];
end


% Salinity limits (plot)
Smin = 31; 
Smax = 37;


% ==========================
[A] = rd_smos_L4_BEC_v1r2(itime,ibasin,data_type);
lon  = A.lon; 
lat = A.lat; 
sss = A.sss;
time_out = A.time_out;

clear A; clc;


lon_min = -180;
lon_max = 180;
lat_min = 40;
lat_max = 90;


[IND] = find(lon >= lon_min & lon <= lon_max & ...
    lat >= lat_min & lat <= lat_max);

lon2 = lon(IND);
lat2 = lat(IND);
sss2 = sss(IND);

xLON = lon_min+0.25:0.25:lon_max-0.25;
yLAT = lat_min-0.25:0.25:lat_max-0.25;

[xLON2,yLAT2] = meshgrid(xLON,yLAT);

sss3 = griddata(lon2,lat2,sss2,xLON2,yLAT2);
%% 
% sss3 = interp2(lon2(iNAN),lat2(iNAN),sss2(iNAN),xLON2,yLAT2,'spline');

lon3 = xLON2;
lat3 = yLAT2;

% # =======================================
%% 01) Plot the contentent of the Arctic dataset (Mercator projection)

% % Salinity limits (plot)
% Smin = 25; 
% Smax = 38;
if plot_type == 1
% [fig 1]
figure; clf
pcolor(lon3,lat3,sss3); shading flat
colormap(jet)
caxis ([Smin Smax])
xlim([xmin xmax]); ylim([-90 90])
fillmap
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(rivers, 'Color', 'blue')
colorbar

title({[basin_str ' SSS ' datestr(itime,'yyyymmdd')];[data_type_str ]})

end
% ========================================
%% 
%% 03) Polar plot projection (follow matlab map limit properties 
%% and axesm documentation)
% 
% Salinity limits (plot)
% Smin = 25; 
% Smax = 38;
if plot_type == 2
load coastlines

% [fig 2]
figure; clf; hold on
axesm('eqaazim','MapLatLimit',[50 90])
axis off
framem on
gridm on
mlabel on
plabel on;
setm(gca,'MLabelParallel',40)

pcolorm(lat3,lon3,sss3); shading flat
colormap(jet)
caxis ([Smin Smax])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1].*0.5)
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(rivers, 'Color', 'blue')
gridm on
title({[basin_str ' SSS ' datestr(itime,'yyyymmdd')];[data_type_str]})

end


%% 03) Lambert plot projection
% Selection Arctic regions of interest (as seen in DUM?)
% 
% 1* Subpolar North Atlantic (40?N 65?N, 90?W 30?E)
% 2* Polar North Atlantic
% 3* North West Arctic Pacific
% 4* North East Arctic Pacific

% 1* Subpolar North Atlantic (40?N 65?N, 90?W 30?E)

if plot_type == 3 && strcmpi(region,'NA') || strcmpi(region,'ALL')

% [fig 3]
figure; clf
lat_min2 = 50;
lat_max2 = 85;
lon_min2 = -100;
lon_max2 = 100;

lat_step = 15;
lon_step = 20;

 Smin2 = 25;
 Smax2 = 38;

map_projection = 'lambert';

fillmap_super(map_projection,lon_min2,lon_max2,...
    lat_min2,lat_max2,lon_step,lat_step)

pcolorm(lat3,lon3,sss3); shading flat
colormap(jet)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(rivers, 'Color', 'blue')
gridm on
title({[basin_str ' North Atlantic (' num2str(lat_min2) ...
    '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
    [data_type_str]})

end

%% 2* Polar North Atlantic (60N 90N, 0 180E)

if plot_type == 3 && strcmpi(region,'NA') || strcmpi(region,'ALL')

% [fig 4]
figure; clf
lat_min2 = 60;
lat_max2 = 85;
lon_min2 = 0;
lon_max2 = 180;

lat_step = 15;
lon_step = 20;

 Smin2 = 33;
 Smax2 = 36;

map_projection = 'lambert';

fillmap_super(map_projection,lon_min2,lon_max2,...
    lat_min2,lat_max2,lon_step,lat_step)


pcolorm(lat3,lon3,sss3); shading flat
colormap(jet)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(rivers, 'Color', 'blue')
gridm on
title({[basin_str ' North Atlantic (' num2str(lat_min2) ...
    '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
    [data_type_str]})

end

%% 3* North West Arctic Pacific


if plot_type == 3 && strcmpi(region,'PA') || strcmpi(region,'ALL')
% [fig 5]
figure; clf
lat_min2 = 45;
lat_max2 = 80;
lon_min2 = 130;
lon_max2 = -120;

lat_step = 15;
lon_step = 20;

Smin2 = 30;
Smax2 = 34;

map_projection = 'lambert';


fillmap_super(map_projection,lon_min2,lon_max2,...
    lat_min2,lat_max2,lon_step,lat_step)

setm(gca, 'Origin', [0 -180]) % Set the Pacific in the middle of the plot

pcolorm(lat,lon,sss); shading flat
colormap(jet)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(rivers, 'Color', 'blue')
gridm on
title({[basin_str ' North Pacific (' num2str(lat_min2) ...
    '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
    [data_type_str]})

end

