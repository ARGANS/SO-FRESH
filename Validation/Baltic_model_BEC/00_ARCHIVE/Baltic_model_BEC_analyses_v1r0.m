% Syntax: Baltic_model_BEC_analyses.m (script)
%
% Description
% Analyses of NEMO (model) data as provided by BEC. The model is:
% "CMEMS V4 Reanalysis: NEMO model 3D fields (monthly means)"
%
%
% Input
%
% For the validation purposes of the Baltic+ there are two periods worth of
% data contained in the following two files:
% *[1] dataset-reanalysis-nemo-monthlymeans_2011_2013.nc
% *[2] dataset-reanalysis-nemo-monthlymeans_2014_2016.nc
%
% The NetCDF files contain the following variables;
%
% [1] latitude and longitude grid geolocation
% [2] mlotst: Mixed Layer Thickness defined by density (as in de Boyer
% Montegut, 2004). -                  MLD   [381x523x36]
% [3] SO: sea_water_salintiy -        SSS   [381x523x36]
% [4] thetao: potential temperature - SST   [381x523x36]
% [5] depth (m)
% [6] time: days since 1950-01-01 00:00:00  [36x1] *to convert to MATLAB number
%
% Baltic Study [2011-2013], contains:
% *[1] Produce Surface Map: SSS, SST and MLD
% *[2] Boxplots ALL-BALTIC, and regional studies
% *[3]
%
% Data source:
% BEC products for validation
%
%
% current version: v1r0 (2020/01/20)
%
% History
% -
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================
clc; clear;
close all


year1 = 2011;    % to choose 2011 or 2014
year2 = year1+2;

year1_str = num2str(year1);
year2_str = num2str(year2);

ibasin = 9;

[xmin,xmax,ymin,ymax,basin_str] = map_lim_raf (ibasin);

% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance from reference point (in Km)


fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/model/']);


save_fig = 1; % flag to save figure [01]; or not [00].
fg_format = 'png';

% create a storage folder
folder_figs = ['/Volumes/Rogue/scratch/Validation/' basin_str '/NEMO/'];

if save_fig == 1
    foldercheck_raf(folder_figs); %! make folder_out
    cd (folder_figs)
    %eval(['!open ' pwd])
end


% ship routes
SHIPS = {'BalticQueen','FinnMaid','Romantica',...
    'SiljaSernade','Transpaper','Victoria'};

% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'NEMO_MISSING_20200121.txt'];


fn = ['dataset-reanalysis-nemo-monthlymeans_' year1_str '_' year2_str];
folder_in = [folder_data];
fn_in = [folder_in fn '.nc'];


% =============
% save vars from netCDF file into the matFILE
param_load = {'latitude','longitude','mlotst','so','thetao','time','depth'};


% NetCDF null constants (for more info see ncdisp(fn)) -- Set by BEC
missing_val = -999;
FILL = missing_val;

% read each variable in the NetCDF file
if exist(fn_in,'file') == 2
    for x = 1:length(param_load)
        param = param_load{x};
        if strcmpi(param,'so')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['SSS_model = squeeze(TT);']);
        elseif strcmpi(param,'thetao')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['SST_model = TT;']);
        elseif strcmpi(param,'LONGITUDE')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['lon = TT;']);
        elseif strcmpi(param,'LATITUDE')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['lat = TT;']);
        elseif strcmpi(param,'mlotst')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['MLD_model = TT;']);
            
        elseif strcmpi(param,'time')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['time_model = TT;']);
        elseif strcmpi(param,'depth')
            TT = double(ncread(fn_in,param));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['depth_model = TT;']);
        end
        clear TT
    end
end

% To meshgrid lon/lat
[lat,lon] = meshgrid(lat,lon);

% Convert the time_model into meaning full reading time
time_num = datenum(time_model);

% time_model :units = "days since 1950-01-01 00:00:00";
time0 = datenum(1950,01,01,0,0,0);
time_out = zeros(size(time_num));

% matlab time number
for n = 1: length(time_num)
    time_out(n) = addtodate(time0,time_num(n),'day');
end

TIME_model = time_out; clear time_out time_model;


%% Plot monthly MAPS
lat_min = ymin;
lat_max = ymax;
lon_min = xmin;
lon_max = xmax;
% Make surface plots of all variables in the model (iparam_str_all) and for
% all TIME_model available
param_str = {'SSS','SST','MLD'};
param_units = {'psu','\circ','m'};


plot_monthly_maps = 1;

if plot_monthly_maps == 1
    for pp = 1:length(iparam_str_all) % loop through variables
        for tt = 1:length(TIME_model)  % loop through time
            %
            iparam_str = param_str{pp};
            itime_num = TIME_model(tt);
            
            eval(['param_in = squeeze(' iparam_str '_model(:,:,' num2str(tt) '));'])
            
            % funcion-snipet to plot baltic model datas
            figure
            [h1,fg_name] = plot_map_model_BEC(lon,lat,param_in,iparam_str,ibasin,itime_num);
            
            title({[iparam_str '_{NEMO} ' basin_str ' ('...
                num2str(lat_min) '-' ...
                num2str(lat_max) '\circ N) ' ]; ...
                [datestr(itime_num,'yyyymmdd')]});
            
            % Save figure  - output -
            folder_this = [folder_figs 'MAPS/' iparam_str '/'];
            
            if fg_save == 1
                foldercheck_raf(folder_this); %! make folder_figs
            end
            
            fg_name = [folder_this fg_name '.' fg_format];
            
            % check fn existence
            fg_exist = exist(fg_name,'file');
            if fg_save == 1 && fg_exist == 0
                save_raf(gcf,fg_name,fg_format); close
            end
        end
    end
    
end


%% Model Regional study (variability within each study region)
% There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
% [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
% [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
% [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
% [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)

% lat and lon as taken from Google Earth

lat_sreg = [...
    55 31 30.22;...
    61 55 50.19;...
    59 35 55.07;...
    57 36 10.05];

lon_sreg = [...
    16 16 15.06;...
    19 14 46.24;...
    23 07 27.96;...
    19 48 25.78];


% convert angles [deg-min-sec] to [deg]
lat_sreg = dms2degrees(lat_sreg);
lon_sreg = dms2degrees(lon_sreg);

SREG = [lat_sreg lon_sreg]; % geolocation centre study regions

SREG_str = {'[01] ArB','[02] BOS','[03] GOF','[04] NBP'};
SREG_str2 = {'ArB','BOS','GOF','NBP'}; % use this for figure names

% Get NEMO-data within a given distance from the center if SREG.
r = 200; % km


%% Plot with all the study regions in the Baltic as seen in DUM (p.28)
figure
plot_study_region(lon_sreg,lat_sreg,r); hold

for nn = 1:length(lon_sreg)
    text(lon_sreg(nn),lat_sreg(nn),sprintf('%02.0f',nn),'fontsize',18)
    
end

xlim([xmin xmax]); ylim([ymin ymax])
box on
grid on

fillmap

title({['Baltic+ Study Regions ']...
    ['(Radius: ' num2str(r) ' km)']});

fg_name = [basin_str '_STUDYREGIONS_NEMO'];

% Save figure  - output -
folder_this = [folder_figs 'MAPS/STUDYREGIONS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name '.' fg_format];

% check fn existence
fg_exist = exist(fg_name,'file');
if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end


%% Get NEMO-data within a searching distance at each study region
[a,b,c] = size(SSS_model);

% number of study region
Nsreg = length(lon_sreg);

% Pre-locate variables
SSS_model_draft = reshape(SSS_model,a*b,c);
SST_model_draft = reshape(SST_model,a*b,c);
MLD_model_draft = reshape(MLD_model,a*b,c);

SSS_irange = nan([a*b,c, Nsreg]);
SST_irange = SSS_irange;
MLD_irange = SSS_irange;

for nn = 1:length(lon_sreg)
    DIST = Distance(lon_sreg(nn),lat_sreg(nn),lon,lat);
    
    irange = find(DIST <= r);
    Ln = length(irange);
    
    SSS_irange(1:Ln,:,nn) = SSS_model_draft(irange,:);
    SST_irange(1:Ln,:,nn) = SST_model_draft(irange,:);
    MLD_irange(1:Ln,:,nn) = MLD_model_draft(irange,:);
    
end


% Return to original shape NEMO-inrange variables
% SSS_irange = reshape(SSS_irange,a,b,c,Nsreg);
% SST_irange = reshape(SST_irange,a,b,c,Nsreg);
% MLD_irange = reshape(MLD_irange,a,b,c,Nsreg);


%% Time series at each location

% 1 > Stats Computation: mean, median, std

% 1.1 >> Baltic Regional stats
SSS_mean_sreg = squeeze(nanmean(SSS_irange,1));
SST_mean_sreg = squeeze(nanmean(SST_irange,1));
MLD_mean_sreg = squeeze(nanmean(MLD_irange,1));

SSS_median_sreg = squeeze(nanmedian(SSS_irange,1));
SST_median_sreg = squeeze(nanmedian(SST_irange,1));
MLD_median_sreg = squeeze(nanmedian(MLD_irange,1));

SSS_std_sreg = squeeze(nanstd(SSS_irange,1));
SST_std_sreg = squeeze(nanstd(SST_irange,1));
MLD_std_sreg = squeeze(nanstd(MLD_irange,1));

% 1.2 >> Baltic-Global stats
SSS_mean_sreg_ALL = squeeze(nanmean(SSS_mean_sreg,2));
SST_mean_sreg_ALL = squeeze(nanmean(SST_mean_sreg,2));
MLD_mean_sreg_ALL = squeeze(nanmean(MLD_mean_sreg,2));

SSS_median_sreg_ALL = squeeze(nanmedian(SSS_median_sreg,2));
SST_median_sreg_ALL = squeeze(nanmedian(SST_median_sreg,2));
MLD_median_sreg_ALL = squeeze(nanmedian(MLD_median_sreg,2));

SSS_std_sreg_ALL = squeeze(nanmean(SSS_std_sreg,2));
SST_std_sreg_ALL = squeeze(nanmean(SST_std_sreg,2));
MLD_std_sreg_ALL = squeeze(nanmean(MLD_std_sreg,2));


%% 2 > Plot stats

% value from 0 (trasparent) to 1 (opaque) (type help boundedline)
shade_transparent = 0.30;

% Loop plots through all the parameteres: SSS, SST and MLD
param_str = {'SSS','SST','MLD'};
param_units = {'psu','\circ','m'};

time_str1 = datestr(TIME_model(1),'yyyymmdd');
time_str2 = datestr(TIME_model(end),'yyyymmdd');


%% 2.1 >> Time Series Mean
figure(22); clf; hold on
for nn = 1 : length(param_str)
    subplot(3,1,nn); hold on
    
    eval(['h1 = plot(TIME_model,'...
        param_str{nn} '_mean_sreg_ALL,'...
        ['''--k'''] ',' ['''linewidth'''] ',2);']); hold on
    eval(['h2 = boundedline(TIME_model,'...
        param_str{nn} '_mean_sreg,'...
        param_str{nn} '_std_sreg_ALL,'...
        ['''alpha'''] ',' ['''orientation'''] ',' ['''vert'''] ','...
        ['''transparency'''] ',' num2str(shade_transparent) ');']); hold on
    
    set(h2,'linewidth',2);
    
    lg = legend([h1;h2],['[00] ALL', SREG_str],'fontsize',14);
    set(lg,'Location','northwest')
    
    box on
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Regional averaged (mean) Baltic ' param_str{nn}],'fontsize',24);
    xlabel('time (months)','fontsize',18);
    ylabel([param_str{nn} ' (' param_units{nn} ')' ], 'fontsize',18);
    
end

fg_name = [basin_str '_MEAN_NEMO_' time_str1 '_' time_str2];

% Save figure  - output -
folder_this = [folder_figs 'TIMESERIES/MEAN/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name '.' fg_format];

% check fn existence
fg_exist = exist(fg_name,'file');
if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end


%% 2.1 >> Time Series Median [box plot at each study region]
for pp = 1: length(SREG_str)
    fg_num = pp*2;
    figure(fg_num); clf; hold on

    
    for nn = 1 : length(param_str)
        subplot(3,1,nn); hold on
        
        eval([ 'boxplot(' param_str{nn} '_irange(:,:,' num2str(pp) '));']); hold on
        
        eval([ 'h1 = plot(' param_str{nn} '_median_sreg(:,' num2str(pp) '),'...
            ['''ok'''] ',' ['''linewidth'''] ',2);'])
        eval([ 'h2 = plot(' param_str{nn} '_mean_sreg(:,' num2str(pp) '),'...
            ['''--k'''] ',' ['''linewidth'''] ',2);'])
        
        lg = legend([h1(1);h2(1)],{' median'; ' mean'},'fontsize',14);
        set(lg,'Location','southwest')

        
        box on
        % datetick('x','mmm-yy','keeplimits','keepticks')
        ylabel([param_str{nn} ' (' param_units{nn} ')' ], 'fontsize',18);
        
        if nn == 1
            title(['Regional averaged (median) Baltic ' SREG_str{pp}],'fontsize',24);
            
        end
        if nn == 3
            xlabel('time (months)','fontsize',18);
        end        
    end
    
    fg_name = [basin_str '_MEADIAN_NEMO_' SREG_str2{pp}...
        time_str1 '_' time_str2];
    
    % Save figure  - output -
    folder_this = [folder_figs 'TIMESERIES/REGIONAL/' SREG_str2{pp} '/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name '.' fg_format];
    
    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
end


%% BOXPLOTS
%% [1] BOXPLOT ALL-BALTIC
iparam_str_all = {'SSS','SST','MLD'};

figure; clf; hold on;
set(gcf,'DefaultAxesFontSize',24);
set(gca,'TickLabelInterpreter','tex');

bplot = boxplot...
    ([SSS_model(:)]);

ylim([0 10]);

ylabel('S')


h = findobj(gca, 'type', 'text');
set(h,'FontSize',25);
set(h,'VerticalAlignment', 'middle');
set(gca,'FontSize',25);

title({['Satellite SSS around Argo float (plat ID:'...
    num2str(platform(nn))...
    ' , PROF#' num2str(ID(nn)) ')'];...
    ['Distance: ' num2str(r) ' km']});

fg_name = [fn '_SAT-SSS_BOXPLOT_PLAT'...
    num2str(platform(nn))...
    '_PROF#' num2str(ID(nn))...
    '_R' num2str(r) 'KM.'...
    fg_format];

% Save figure  - output -
folder_this = [folder_figs 'BOXPLOTS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];


fg_exist = exist(fg_name,'file'); % check fn existence

if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end





