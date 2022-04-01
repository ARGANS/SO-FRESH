% Plot SSS Baltic
% 
% Description
% Plot monthly climatology from SSS (BEC) data in the Baltic
% 
% version: v1r1
% This version uses the BEC Global SSS data. 
% Future releases contemplate
% using the regional Baltic+ Salinity data (to be ready in September 2019) 
%
% 
% =========================================================================
% Date
% 13/08/2019 
%
% rcatany
%
% 
% =========================================================================
%%
% clc;
% clear

% script switches
% DATA_PATH = '/Users/rejc1e11/Documents/Data/';    % # Use this path working local
DATA_PATH = '/Volumes/Rogue/Data/';                 % # Use this path working hard drive
DATA_PRODUCT = 'Global';                            % # BEC product type: [1] Global; [2] Arctic (avail. Sep 2019); [3] Baltic (avail Sep 2019)
folder_fig = '/Users/rejc1e11/Documents/WORK/Argans/PROJECTS/Pictures/';
data_type = 3;

plot_example = 1;   % plot example
save_fig = 1;       % Save example figures
fig_format = 'png';

% Salinity limits (plot)
Smin = 0; 
Smax = 15;

% # BEC product type: [1] Global; [2] Arctic (avail. Sep 2019); 
% [3] Baltic (avail Sep 2019)

% ARC_PRODUCT: flag (conditional) to know when Arctic regional product is in
% use; Values of ARC_PRODUCT: [1] use arctic product; [0] use global
% product
product_global_type = [3,4]; %set flag to identify when to use the global product
GLO_PRODUCT = ibasin ~= 7 || ismember(data_type,product_global_type);

if data_type == 1 % Arctic+ product (v3.0)
    data_type_str = '9_days';
    data_version_str = 'Arctic_v3.0';
    ndays = 9;
    
    ARC_PRODUCT = ibasin == 7 &&  data_type == 1;% Arctic+ product (v3.0)
    folder_data = 'SSS/Arctic/BEC/v3/30_beta_12/'; 
    
elseif data_type == 2 % Arctic BEC product (v2.0)
    data_type_str = '9_days';
    data_version_str = 'Arctic_v2.0';
    ndays = 9;
    
    ARC_PRODUCT = ibasin == 7 && data_type == 2; % Arctic BEC prodcut (v2.0)
    folder_data = 'SSS/Arctic/BEC/v2/';

elseif data_type == 3 % Use the BEC Global product (v001)
    data_type_str = 'global';
    data_version_str = 'Global_v001';
    ndays = 9;
    
    ARC_PRODUCT = 0;
    folder_data = 'SSS/Baltic/BEC/GlobalProduct/';

else % Use the BEC Global product (default mode)
    data_type_str = 'global';
    data_version_str = 'Global_v001';
    ndays = 9;
    ARC_PRODUCT = 0;
    folder_data = 'SSS/Global/BEC/L4/';
end


% Example (settings)
ibasin = 9; % Basin number (Baltic, ibasin 9)

iyear = 2011;
imonth = 2;
iday = 15;
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

lon_min = xmin;
lon_max = xmax;
lat_min = ymin;
lat_max = ymax;


% =================================
% Load monthly climatology (2011-2017)
folder_anom = [DATA_PATH folder_data 'Climatology/'];
fn = [folder_anom 'Baltic_CLIMAT_20112016_' sprintf('%02.0f',imonth) '.mat'];
load(fn);

SSS_CLIMAT = mean(SSS_climat_mean,3);


% =================================
% Load monthly data
folder_monthly = [DATA_PATH folder_data 'monthly/'...
    sprintf('%02.0f',iyear) '/'];
fn = [folder_monthly basin_str '_MONTHLY_' sprintf('%02.0f',iyear)...
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


% # =======================================
%% 01) Plot the contentent of the Baltic dataset (Mercator projection)
lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;


figure (1); clf
pcolor(lon3,lat3,sss3); shading flat
colormap(jet);
caxis ([Smin Smax])
xlim([lon_min2 lon_max2]); ylim([lat_min2 lat_max2])
fillmap([0.5 0.5 0.5],1)
cbar = colorbar;

title([basin_str ' North Atlantic (Mercator Projection ' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'yyyymmdd')])
xlabel('longitude [?]');
ylabel('latitude [?]');
ylabel (cbar,'S (psu)');


if plot_example == 1 && save_fig == 1
    folder_fig_out =[folder_fig 'Projections/'];
    foldercheck_raf(folder_fig_out)
    
    fig_name = [folder_fig_out basin_str '_MERCATOR_' datestr(iTIME,'yyyymmdd') '.' fig_format];
    fig_exist = exist(fig_name,'file');
    
    if save_fig == 1 && fig_exist == 0
        save_raf(gcf,fig_name,fig_format)
    end
end



% ========================================
%% 
%% 02) Polar plot projection (follow matlab map limit properties and axesm documentation)
% 
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

if plot_example == 1 && save_fig == 1
    folder_fig_out =[folder_fig 'Projections/'];
    foldercheck_raf(folder_fig_out)
    
    fig_name = [folder_fig_out basin_str '_POLAR_' datestr(iTIME,'yyyymmdd') '.' fig_format];
    fig_exist = exist(fig_name,'file');
    
    if save_fig == 1 && fig_exist == 0
        save_raf(gcf,fig_name,fig_format)
    end
end

%% 03) Lambert plot projection
% Selection Arctic regions of interest (as seen in DUM?)
%% 
% * Subpolar North Atlantic (40?N 65?N, 90?W 30?E)
% * Polar North Atlantic (65?N 90?N, 40?W 50?E) 
% * North West Arctic Pacific
% * North East Arctic Pacific
% Subpolar North Atlantic (40?N 65?N, 90?W 30?E) [Notice line OSNAP array]
% 


figure(3); clf
lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;

lat_step = 5;
lon_step = 10;



axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
    'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
    'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)

pcolorm(lat3,lon3,sss3); shading flat
colormap(jet)
caxis ([Smin Smax])
colorbar
framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title([basin_str ' North Atlantic (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'yyyymmdd')])
hold on


if plot_example == 1 && save_fig == 1
    folder_fig_out =[folder_fig 'Projections/'];
    foldercheck_raf(folder_fig_out)
    
    fig_name = [folder_fig_out basin_str '_LAMBERT_' datestr(iTIME,'yyyymmdd') '.' fig_format];
    fig_exist = exist(fig_name,'file');
    
    if save_fig == 1 && fig_exist == 0
        save_raf(gcf,fig_name,fig_format)
    end
end

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
lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;

lat_step = 5;
lon_step = 10;

Smin2 = 0;
Smax2 = 15;

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
%% Baltic SSS Monthly Climatology (2011-2017)
figure(4); clf
lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;

lat_step = 15;
lon_step = 20;

Smin2 = 0;
Smax2 = 15;

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



if plot_example == 1 && save_fig == 1
    folder_fig_out =[folder_fig 'CLIMAT/'];
    foldercheck_raf(folder_fig_out)
    
    fig_name = [folder_fig_out basin_str 'CLIMAT_MONTHLY_' datestr(iTIME,'yyyymmdd') '.' fig_format];
    fig_exist = exist(fig_name,'file');
    
    if save_fig == 1 && fig_exist == 0
        save_raf(gcf,fig_name,fig_format); close
    end
end


%% Baltic SSS (daily) Anomaly
SSS_delta_daily = sss3 - SSS_CLIMAT3;

figure(5); clf
lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;

lat_step = 5;
lon_step = 5;

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
title([basin_str ' (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'yyyymmdd')])
hold on




if plot_example == 1 && save_fig == 1
    folder_fig_out =[folder_fig 'ANOMALIES/'];
    foldercheck_raf(folder_fig_out)
    
    fig_name = [folder_fig_out basin_str 'ANOM_DAILY_' datestr(iTIME,'yyyymmdd') '.' fig_format];
    fig_exist = exist(fig_name,'file');
    
    if save_fig == 1 && fig_exist == 0
        save_raf(gcf,fig_name,fig_format); close
    end
end


%% Baltic SSS (monthly) Anomaly

SSS_delta_monthly = sss_monthly3 - SSS_CLIMAT3;

figure(5); clf
lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;

lat_step = 5;
lon_step = 5;

Smin2 = -1;
Smax2 = 1;

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
title([basin_str ' (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) SSS ' datestr(iTIME,'mmm-yyyy')])
hold on



if plot_example == 1 && save_fig == 1
    folder_fig_out =[folder_fig 'ANOMALIES/'];
    foldercheck_raf(folder_fig_out)
    
    fig_name = [folder_fig_out basin_str 'ANOM_MONTHLY_' datestr(iTIME,'yyyymmdd') '.' fig_format];
    fig_exist = exist(fig_name,'file');
    
    if save_fig == 1 && fig_exist == 0
        save_raf(gcf,fig_name,fig_format); close
    end
end




%% 
%