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
% [1] dataset-reanalysis-nemo-monthlymeans_2011_2013.nc
% [2] dataset-reanalysis-nemo-monthlymeans_2014_2016.nc
%
% Based on the availability of data, the inputs of the scripts are:
% *[0] time_range = 2011_2016 * ALL data available
% *[1] time_range = 2011_2013 * Use first time range
% *[2] time_range = 2014_2016 * Use second time range
%
%
% The NetCDF files contain the following variables;
%
% [1] latitude and longitude grid geolocation
% [2] mlotst: Mixed Layer Thickness defined by density (as in de Boyer
% Montegut, 2004). -                  MLD   [381x523x36]
% [3] SO: sea_water_salintiy -        SSS   [381x523x36]
% [4] thetao: potential temperature - SST   [381x523x36]
% [5] depth (m)
% [6] time: days since 1950-01-01 00:00:00  [36x1]*to convert to MATLAB number
%
% Baltic Study, contains: *limited by SSS-SMOSdata
% *[1] Produce Surface Map: SSS, SST and MLD        [2011-2016]
% *[2] Boxplots ALL-BALTIC, and regional studies    [2011-2016]
% *[3] Comparison SSS-NEMO Vs SSS-SMOS              [2011-2013] 
%
% Data source:
% BEC products for validation
%
%
% current version: v1r1 (2020/01/29)
%
% History
% v1r0 | 20200120 - script creation
% v1r1 | 20200119 - use nemo_rd_BEC_Baltic_FUN.m , and shorten this script
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================
clc; clear;
close all

% Inputs-script
time_range = 0; % [0] 2011-2016; [1] 2011-2013; [2] 2014-2016
ibasin = 9;     % [9] Baltic


%% snippet to load model_data
% this snipped uses function 
% [data_model] = nemo_rd_BEC_Baltic_FUN (iyear,imonth,iday,ibasin,grid_output)
grid_output = 1; % grid nemo to SMOS-BEC grid [1]; or not [0]; 

Baltic_model_BEC_load_data; % snipped


% Get Baltic limits
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

%% Plot monthly MAPS
lat_min = ymin;
lat_max = ymax;
lon_min = xmin;
lon_max = xmax;
% Make surface plots of all variables in the model (iparam_str_all) and for
% all TIME_model available
param_str = {'SSS','SST','MLD'};
param_units = {'psu','\circ','m'};

plot_monthly_maps = 0;

if plot_monthly_maps == 1
    for pp = 1:length(param_str) % loop through variables
        for tt = 1:length(TIME_model)  % loop through time
            %
            iparam_str = param_str{pp};
            itime_num = TIME_model(tt);
            
            eval(['param_in = squeeze(' iparam_str...
                '_model(:,:,' num2str(tt) '));'])
            
            % funcion-snipet to plot baltic model datas
            figure
            [h1,fg_name] = plot_map_model_BEC(lon,lat,...
                param_in,iparam_str,ibasin,itime_num);
            
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
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(rivers, 'Color', 'blue')

xlabel('longitude (\circ)','fontsize',24)
ylabel('latitude (\circ)','fontsize',24)
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


%% Time series at each location

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
    ylabel([param_str{nn} ' (' param_units{nn} ')' ], 'fontsize',24);

    % Xticklabels only in the bottom plot of the figure
    if nn == 1 || nn == 2
        set(gca,'xticklabel',{[]})
    elseif nn == 3
        datetick('x','mmm-yy','keeplimits','keepticks')
        xlabel('time (months)','fontsize',24);
        
    end
    
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

N_int_labels = 1; % time interval (number of months)

tick_numbers = 1:N_int_labels:length(TIME_model);

tick_labels = {tick_numbers};

for pp = 1 : length(SREG_str) % loop each study region
    fg_num = pp*2;
    figure(fg_num); clf; hold on
    
    
    for nn = 1 : length(param_str) % loop each parameter
        subplot(3,1,nn); hold on
        
        eval(['boxplot(' param_str{nn} '_irange(:,:,' num2str(pp) '));']); hold on
        
        eval(['h1 = plot(' param_str{nn} '_median_sreg(:,' num2str(pp) '),'...
            ['''ok'''] ',' ['''linewidth'''] ',2);'])
        eval(['h2 = plot(' param_str{nn} '_mean_sreg(:,' num2str(pp) '),'...
            ['''--k'''] ',' ['''linewidth'''] ',2);'])
        
        lg = legend([h1(1);h2(1)],{' median'; ' mean'},'fontsize',14);
        set(lg,'Location','southwest')
        
        box on
        ylabel([param_str{nn} ' (' param_units{nn} ')' ], 'fontsize',18);
        
        
        xticks(tick_numbers)
        xticklabels(tick_labels{:})

        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
            
        if nn == 1
            title(['Regional averaged (median) Baltic '...
                SREG_str{pp}],'fontsize',24);
        end
        
        
        if nn == 1 || nn == 2
            set(gca,'xticklabel',{[]})
        elseif nn == 3
            % datetick('x','mmm-yy','keeplimits','keepticks')
            xticks(tick_numbers)
            xticklabels(tick_labels{:})
            xlabel('time (months)','fontsize',24);
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


%% [3] SSS-NEMO Vs SSS-SMOS (Baltic+ data)

%% [3.1] Get SSS-SMOS data (Baltic+ data)

% % *[1] Baltic BEC NOM product (v1.0)
% data_type = 4;  % [4] Baltic+ NOM (v1.0), [5] Baltic NS (v1.0), [6] Global BEC
% plot_type = 3;  % [1] Mercator, [2] Polar, [3] Lambert (conical), [4] All-plots
% region = 'Baltic';
% Baltic_plot_fun(ibasin,itime,data_type,plot_type,region)
% 
% % *[2] Baltic+ product (v3.0)
% data_type = 5;  % [4] Baltic+ NOM (v1.0), [5] Baltic NS (v1.0), [6] Global BEC
% plot_type = 3;  % [1] Mercator, [2] Polar, [3] Lambert (conical), [4] All-plots
% region = 'Baltic';
% Baltic_plot_fun(ibasin,itime,data_type,plot_type,region)


