% Plot SSS Arctic
% 
% Description
% Plot to see SSS data in the Arctic
% v1r1
%
%
% Date
% 31/03/2019 
%
% rcatany
%
% 
% =======================================================================
%%

% clc;
% clear

% script switches
% DATA_PATH = '/Users/rejc1e11/Documents/Data/'; % # Use this path working local
DATA_PATH = '/Volumes/Rogue/Data/';            % # Use this path working hard drive

folder_fig = ['/Users/rejc1e11/Documents/WORK/Argans/PROJECTS/Pictures/'];


plot_example = 1;
save_fig = 1;
fig_format = 'png';


% Example (settings)
ibasin = 9; % Basin number (Baltic, ibasin 9)

iyear = 2016;
imonth = 7;
iday = 20;
iTIME = datenum(iyear,imonth,iday);
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

folder_fig = [folder_fig basin_str '/'];
foldercheck_raf(folder_fig)



% ==========================
[A] = rd_smos_L4_BEC_v1r1(datenum(iyear,imonth,iday),ibasin);
lon  = A.lon; 
lat = A.lat; 
sss = A.sss;
time_out = A.time_out;

clear A; clc;


lon_min = -180;
lon_max = 180;
lat_min = 40;
lat_max = 90;


% =================================
% Load monthly climatology (2011-2017)
folder_anom = [DATA_PATH 'SSS/Arctic/BEC/Climatology/'];
fn = [folder_anom 'Arct_CLIMAT_20112017_' sprintf('%02.0f',imonth) '.mat'];
load(fn);

SSS_CLIMAT = mean(SSS_climat_mean,3);


% =================================
% Load monthly data
folder_monthly = [DATA_PATH 'SSS/Arctic/BEC/monthly/' sprintf('%02.0f',iyear) '/'];
fn = [folder_monthly 'Arctic_MONTHLY_' sprintf('%02.0f',iyear)...
    '_' sprintf('%02.0f',imonth) '.mat'];
load(fn)
sss_monthly = sss_mean;


% ================================
% Griddata to homogeneous grid

[IND] = find(lon >= lon_min & lon <= lon_max & lat >= lat_min & lat <= lat_max);

lon2 = lon(IND);
lat2 = lat(IND);
sss2 = sss(IND);
SSS_CLIMAT2 = SSS_CLIMAT(IND);
sss_monthly2 = sss_monthly(IND);


xLON = lon_min+0.25:0.25:lon_max-0.25;
yLAT = lat_min-0.25:0.25:lat_max-0.25;

[xLON2,yLAT2] = meshgrid(xLON,yLAT);
lon3 = xLON2;
lat3 = yLAT2;

sss3 = griddata(lon2,lat2,sss2,xLON2,yLAT2);
SSS_CLIMAT3 = griddata(lon2,lat2,SSS_CLIMAT2,xLON2,yLAT2);
sss_monthly3 = griddata(lon2,lat2,sss_monthly2,xLON2,yLAT2);
%% 
% 


% # =======================================
%% 01) Plot the contentent of the Arctic dataset (Mercator projection)

% Salinity limits (plot)
Smin = 25; 
Smax = 38;

lat_min2 = -90;
lat_max2 = 90;
lon_min2 = -180;
lon_max2 = 180;


figure (1); clf
pcolor(lon3,lat3,sss3); shading flat
colormap(jet);
caxis ([Smin Smax])
xlim([lon_min2 lon_max2]); ylim([-90 90])
fillmap
cbar = colorbar;

title([basin_str ' North Atlantic (Mercator Projection ' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'yyyymmdd')])
xlabel('longitude (?)');
ylabel('latitude (?)');
ylabel (cbar,'S (psu)')
 
% folder_fig_out =[folder_fig 'Projections/'];
% foldercheck_raf(folder_fig_out)
%
% fig_name = [folder_fig_out 'Arctic_MERCATOR_' datestr(iTIME,'yyyymmdd')];
% fig_exist = exist(fig_name,'file');
% 
% if save_fig == 1 && fig_exist == 0
%     save_raf(gcf,fig_name,fig_format)
% end




% ========================================
%% 
%% 02) Polar plot projection (follow matlab map limit properties and axesm documentation)
% 
% Salinity limits (plot)
Smin = 25; 
Smax = 38;

load coastlines

figure (2); clf; hold on
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
gridm on
title([basin_str ' SSS ' datestr(iTIME,'yyyymmdd')])



%% 03) Lambert plot projection
% Selection Arctic regions of interest (as seen in DUM?)
%% 
% * Subpolar North Atlantic (40?N 65?N, 90?W 30?E)
% * Polar North Atlantic (65?N 90?N, 40?W 50?E) 
% * North West Arctic Pacific
% * North East Arctic Pacific
% Subpolar North Atlantic (40?N 65?N, 90?W 30?E) [Notice line OSNAP array]
% 
% Get OSNAP location
[OSNAP] = osnap_rd_v1r1('TS')

lon_osnap = OSNAP.lon;
lat_osnap = OSNAP.lat;

figure(3); clf
lat_min2 = 45;
lat_max2 = 80;
lon_min2 = -100;
lon_max2 = 30;

lat_step = 15;
lon_step = 20;

Smin2 = 33;
Smax2 = 36;

axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
    'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
    'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)

pcolorm(lat3,lon3,sss3); shading flat
colormap(jet)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title([basin_str ' North Atlantic (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'yyyymmdd')])
hold on

plotm(lat_osnap,lon_osnap,'ok')



%% Subpolar North Atlantic (40?N 65?N, 90?W 30?E) [Notice OSCAR currents]
% 
fn = [DATA_PATH 'Currents/oscar_vel' num2str(iyear) '.nc'];
[lat_oscar, lon_oscar, time_oscar, date_oscar, depth_oscar, mask, u, v, uf, vf] = read_OSCAR(fn);


time_ind = find(time_oscar >= iTIME,1,"first"); % index the time for current plot (iTIME)

time_oscar = time_oscar(time_ind);
[LON_oscar,LAT_oscar] = meshgrid(lat_oscar,lon_oscar);
U = squeeze(u(:,:,1,time_ind));
V = squeeze(v(:,:,1,time_ind));

%%

figure(4); clf
lat_min2 = 45;
lat_max2 = 80;
lon_min2 = -100;
lon_max2 = 30;

lat_step = 15;
lon_step = 20;

Smin2 = 33;
Smax2 = 36;

axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
    'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
    'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)

pcolorm(lat3,lon3,sss3); shading flat
colormap(jet)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title({[basin_str ' North Atlantic (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'yyyymmdd')],['OSCAR currents ( m/s))']})
hold on

h1 = quiverm(LON_oscar,LAT_oscar,V,U,'k',10);
% h1 = quivermc(LAT_oscar,LON_oscar,U,V,'color','k','reference',10,'linewidth',2);
% set(h1,'AutoScale','on', 'AutoScaleFactor', 2)




%% Computation of SSS anomalies
% % 
% If we want to analyse how unusual a particular event or time period is, it?s 
% helpful to examine the anomaly from the long-term mean rather than the values 
% themselves. 
% 
% The long-term mean is referred to as the climatology. This is the mean of 
% each month with the year range (i.e. the average January, average February, 
% average March etc.).
% 
% The *anomaly is the value minus the respective climatology value*.  So the 
% anomaly of SSS in August 2011 is: August 2011 minus the average of all the Augusts 
% in the time series.
%% 
%% Arctic SSS Monthly Climatology (2011-2017)

figure(4); clf
lat_min2 = 45;
lat_max2 = 80;
lon_min2 = -100;
lon_max2 = 30;

lat_step = 15;
lon_step = 20;

Smin2 = 30;
Smax2 = 38;

axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
    'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
    'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)

pcolorm(lat3,lon3,SSS_CLIMAT3); shading flat
colormap(jet)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title([basin_str 'Monthly climatology (2011-2017) North Atlantic (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'mmm')])
hold on

plotm(lat_osnap,lon_osnap,'ok')



%% Arctic SSS (daily) Anomaly


SSS_delta_daily = sss3 - SSS_CLIMAT3;

figure(5); clf
lat_min2 = 45;
lat_max2 = 80;
lon_min2 = -100;
lon_max2 = 30;

lat_step = 15;
lon_step = 20;

Smin2 = -0.5;
Smax2 = 0.5;

axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
    'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
    'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)

pcolorm(lat3,lon3,SSS_delta_daily); shading flat
colormap(french2)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title([basin_str ' North Atlantic (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) ?SSS ' datestr(iTIME,'yyyymmdd')])
hold on

plotm(lat_osnap,lon_osnap,'ok')




%% Arctic SSS (monthly) Anomaly

SSS_delta_monthly = sss_monthly3 - SSS_CLIMAT3;

figure(5); clf
lat_min2 = 45;
lat_max2 = 80;
lon_min2 = -100;
lon_max2 = 30;

lat_step = 15;
lon_step = 20;

Smin2 = -0.5;
Smax2 = 0.5;

axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
    'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
    'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)

pcolorm(lat3,lon3,SSS_delta_monthly); shading flat
colormap(french2)
caxis ([Smin2 Smax2])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title([basin_str ' North Atlantic (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) ?SSS ' datestr(iTIME,'mmm-yyyy')])
hold on

plotm(lat_osnap,lon_osnap,'ok')





%% 
%