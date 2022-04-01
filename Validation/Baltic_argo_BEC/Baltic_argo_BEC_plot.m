% Baltic_argo_BEC_plot
%
% Description
% Plot Argo floats as provided by BEC.
% Input: matlab files created from the original netCDF files (use
% ar_rd_BEC.m)
%
% Data source:
% BEC products for validation
%
%
%
% v1r0 (2019/10/31)
% History
% [1] creation of this script [20191031]
%
% ========================================================================

clc; clear;

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/argo/argo_mat/']);

iyear = 2013;
imonth = 10;
iday = 1;
% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file

ibasin = 9; % Basin number: 9 (Baltic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);


folder_in = ([folder_data...
    num2str(iyear) '/' sprintf('%02.0f',imonth) '/']);


itime_start = datenum(iyear,imonth,iday);
itime_end = itime_start + ndays;

% load matlab file [fn]
fn = (['argo_' ...
    datestr(itime_start,'yyyymmdd') '_'...
    datestr(itime_end,'yyyymmdd')]);

fn_in = ([folder_in fn '.mat']);

if exist(fn_in,'file') == 2
    load (fn_in)
elseif exist(fn_in,'file') == 0
    error (['file: ' fn_in ' needs to be downloaded'])
end

% number of profiles (observations)
Nobs = length(lon);


%% *[1] Plot Argo geolocation

figure (1); clf; hold on

lat_min2 = ymin;
lat_max2 = ymax;
lon_min2 = xmin;
lon_max2 = xmax;

lat_step = 50;
lon_step = 90;

load coastlines

axesm('eqaazim','MapLatLimit',[50 90])
axis off
framem on
gridm on
mlabel on
plabel on;
setm(gca,'MLabelParallel',40,'Fontsize',14)

framem on
axis off
tightmap
geoshow('landareas.shp','facecolor',[1 1 1]*.5)
gridm on
title({[basin_str ' (' num2str(lat_min2) '-' ...
    num2str(lat_max2) '\circ N) ' ]; ...
    ['N_{BEC} = ' num2str(Nobs)];...
    [datestr(itime_start,'yyyymmdd') '-'...
    datestr(itime_end,'yyyymmdd')]});
hold on

plotm(lat,lon,'ko','MarkerFaceColor','r')


%% *[2] Plot Argo profiles

%% Interpolate Z-direction (narrow depth levels)
% Interpolation specs: Grid dimensions
pres1 = 0;              % min pres bin
pres2 = max(PRES);            % max depth (m)
grid_p = 0.5;           % interpolant distance (m) between depth levels

int_method = 'linear';              % interpolation method (see doc interp1)
pres_p = [pres1:grid_p:pres2]';     % interpolated presure
zlevels = length(pres_p);           % number of depth levels in the interpolated vars

[i,j] = size(TEMP);

t2 = nan(zlevels,j);
s2 = nan(zlevels,j);
p2 = ones(zlevels,j).*pres_p;

PRES2 = repmat(PRES,[1,j]);

% interpolate profiles z-axis (depth)
for j = 1 : Nobs
    t1 = TEMP(:,j);
    s1 = PSAL(:,j);
    p1 = PRES2(:,j); % presure levels are the same for all the profiles
    
    ind = isnan(t1) & isnan(p1) & isnan(s1) | p1 > pres2;
    
    p1(ind) = [];
    t1(ind) = [];
    s1(ind) = [];
    
    [p1, index] = unique(p1);
    
    
    dp = diff(p1);
    
    if min(dp) >= 0.01
        t2_interp = interp1(p1,t1(index),pres_p,int_method);
        t2(:,j) = t2_interp;
        
        s2_interp = interp1(p1,s1(index),pres_p,int_method);
        s2 (:,j) = s2_interp;
    end
    % clear vars within loop
    clear t1 s1 p1 ind
end

t_intp = t2;
s_intp = s2;
p_intp = p2;

clear t2 s2 p2


%% Example
%% 1/ plot raw Vs interp T and S of a given Profile (nprof)
nprof = find(platform==6901150); % Location North Atlantic[20161101]
%nprof = find(platform == 4901821); % Location: East Pacific [20161101];

if ~isempty(nprof)
    
    T = TEMP(:,nprof);
    S = PSAL (:,nprof);
    P = PRES2 (:,nprof);
    
    T_intp = t_intp(:,nprof);
    S_intp = s_intp(:,nprof);
    P_intp = p_intp(:,nprof);
    
    plat = platform(nprof);
    
    Tmin = 0;
    Tmax = 15;
    Smin = 0;
    Smax = 40;
    
    figure
    subplot(1,2,1)
    h1 = plot (T, -P,'k--');
    hold on
    h2 = plot (T_intp, -P_intp,'ko','MarkerFaceColor','r');
    
    ylim([-pres2 -pres1])
    %xlim ([Tmin Tmax])
    
    xlabel ('T (\circ)'); ylabel('depth (m)')
    grid on
    legend ([h1;h2],{'raw';'interp'},'location','SouthEast')
    title(['plaform: ' num2str(plat)]);
    
    subplot(1,2,2)
    h1 = plot (S, -P,'k--');
    hold on
    h2 = plot (S_intp(:,1), -P_intp,'ko','MarkerFaceColor','b');
    
    ylim([-pres2 -pres1])
    %xlim ([Smin Smax])
    
    xlabel ('S (psu)'); ylabel('depth (m)')
    grid on
    legend ([h1;h2],{'raw'  ;'interp'},'location','SouthWest')
    title(['plaform: ' num2str(plat)]);
    
end

%% 2/ Show Geolocation of a given profile (nprof)
figure (1); hold on;

ind1 = lon >=-100 & lon <=30; % floats in NA
ind2 = lon >= 100 & lon <=180; % floats in PA


h1 = plotm(lat(ind1),lon(ind1),'ko','MarkerFaceColor','r');
h2 = plotm(lat(ind1),lon(ind1),'ko','MarkerFaceColor','b');

if ~isempty(nprof)
    
    h3 = plotm(lat(nprof),lon(nprof),'pk','MarkerFaceColor','y','MarkerSize',20);
    
end


