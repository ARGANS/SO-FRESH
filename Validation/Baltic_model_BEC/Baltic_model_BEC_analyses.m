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
time_range = 1;
% *[0] time_range = 2011_2016 * ALL data available
% *[1] time_range = 2011_2013 * Use first time range
% *[2] time_range = 2014_2016 * Use second time range


ibasin = 9;

%% snippet to load model_data
% this snipped uses function
% [data_model] = nemo_rd_BEC_Baltic_FUN (iyear,imonth,iday,ibasin,grid_output)
grid_output = 1; % grid nemo to SMOS-BEC grid [1]; or not [0];

Baltic_model_BEC_load_data;


% Get Baltic limits
[xmin,xmax,ymin,ymax,basin_str] = map_lim_raf (ibasin);

% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance from reference point (in Km)

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/model/']);

fg_save = 1; % flag to save figure [01]; or not [00].
fg_format = 'png';

% create a storage folder
folder_figs = ['/Volumes/Rogue/scratch/Validation/' basin_str '/NEMO/'];

if fg_save == 1
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
param_model_str = {'SSS_{model}','SST_{model}','MLD_{model}'};
param_units = {'psu','\circ','m'};

plot_monthly_maps = 0;

if plot_monthly_maps == 1
    for pp = 1:length(param_model_str) % loop through variables
        for tt = 1:length(TIME_model)  % loop through time
            %
            iparam_str = param_model_str{pp};
            itime_num = TIME_model(tt);
            
            eval(['param_in = squeeze(' iparam_str...
                '(:,:,' num2str(tt) '));'])
            
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
plot_time_series = 0;
%% 2 > Plot stats

% value from 0 (trasparent) to 1 (opaque) (type help boundedline)
shade_transparent = 0.30;

% Loop plots through all the parameteres: SSS, SST and MLD
param_model_str = {'SSS_model','SST_model','MLD_model'};
param_units = {'psu','\circ','m'};

time_str1 = datestr(TIME_model(1),'yyyymmdd');
time_str2 = datestr(TIME_model(end),'yyyymmdd');

% string for labels (sub-indices)
param_model_str2 = {'SSS_{model}','SST_{model}','MLD_{model}'};


if plot_time_series == 1
    %% 2.1 >> Time Series Mean
    figure(22); clf; hold on
    set(gcf,'DefaultAxesFontSize',24);
    for nn = 1 : length(param_model_str)
        subplot(3,1,nn); hold on
        
        eval(['h1 = plot(TIME_model,'...
            param_model_str{nn} '_mean_sreg_ALL,'...
            ['''--k'''] ',' ['''linewidth'''] ',2);']); 
        hold on
        
        eval(['h2 = boundedline(TIME_model,'...
            param_model_str{nn} '_mean_sreg,'...
            param_model_str{nn} '_std_sreg_ALL,'...
            ['''alpha'''] ',' ['''orientation'''] ',' ['''vert'''] ','...
            ['''transparency'''] ','...
            num2str(shade_transparent) ');']); 
        hold on
        
        set(h2,'linewidth',2);
        
        lg = legend([h1;h2],['[00] ALL', SREG_str],'fontsize',14);
        set(lg,'Location','northwest')
        
        grid on
        grid minor
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        
        % Xticklabels only in the bottom plot of the figure
        if nn == 1 || nn == 2
            set(gca,'xticklabel',{[]})
        elseif nn == 3
            datetick('x','mmm-yy','keeplimits','keepticks')
            xlabel('time (months)','fontsize',24);
            
        end
        
        title(['Regional averaged (mean) Baltic '...
            param_model_str2{nn}],'fontsize',24);
        ylabel([param_model_str2{nn}...
            ' (' param_units{nn} ')' ], 'fontsize',24);
        
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
        set(gcf,'DefaultAxesFontSize',24);
        
        for nn = 1 : length(param_model_str) % loop each parameter
            subplot(3,1,nn); hold on
            
            eval(['boxplot(' param_model_str{nn} '_irange(:,:,' num2str(pp) '));']); hold on
            
            eval(['h1 = plot(' param_model_str{nn} '_median_sreg(:,' num2str(pp) '),'...
                ['''ok'''] ',' ['''linewidth'''] ',2);'])
            eval(['h2 = plot(' param_model_str{nn} '_mean_sreg(:,' num2str(pp) '),'...
                ['''--k'''] ',' ['''linewidth'''] ',2);'])
            
            lg = legend([h1(1);h2(1)],{' median'; ' mean'},'fontsize',14);
            set(lg,'Location','southwest')
            
            box on
            ylabel([param_model_str{nn} ' (' param_units{nn} ')' ], 'fontsize',24);
            
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
    
end % plot_time_series





%% [3] SSS-NEMO Vs SSS-SMOS (Baltic+ data)

if time_range == 0
    clc
    disp('To see the SSS-SMOS Vs NEMO, use: Baltic_modelvsSMOS.m (script)')
    disp ('>> NOTE: model vs SMOS only available in time_range [1]: 2011-2013')
    
elseif time_range ~= 0 % Use modelVsSatellite script to plot this section
    
    
    %% [3.1] Get SSS-SMOS data (Baltic+ data)
    
    % Baltic+ [v1.0] data ranges from: 20110201 to 20131227
    iyear1 = 2011;
    imonth1 = 02;
    iday1 = 01;
    
    iyear2 = 2013;
    imonth2 = 12;
    iday2 = 27;
    
    itime_min = datenum(iyear1,imonth1,iday1,0,0,0);
    itime_max = datenum(iyear2,imonth2,iday2,0,0,0);
    
    time_range = itime_min:itime_max;
    
    % Convert to squared matrices lon_irange and lat_irange (NEMO data in each
    % region)
    [a,b,c] = size(SSS_model);
    Nsreg = 4;
    lon_irange = reshape(lon_irange,a,b,Nsreg);
    lat_irange = reshape(lat_irange,a,b,Nsreg);
    
    folder_out = ['/Volumes/Rogue/scratch/Validation/Baltic/NEMO/tmp/'];
    foldercheck_raf(folder_out);
    
    fn_out = [folder_out 'Baltic_SSS_REGIONAL_'...
        datestr(itime_min,'yyyymmdd')...
        '_' datestr(itime_max,'yyyymmdd') '.mat'];
    
    
    % Get SSS-SMOS data in each Baltic-region and to same grid size as
    % model_data
    
    fn_exist = exist(fn_out,'file');
    
    vars_save = {...
        'SSS_NM','SSSe_NM','SSS_NS','SSSe_NS',...
        'SSS_NM_irange','SSSe_NM_irange',...
        'SSS_NS_irange','SSSe_NS_irange',...
        'lon_smos','lat_smos',...
        'lon_sreg','lat_sreg','time_range'};
    
    
    if fn_exist ~= 2
        for tt = 1:length(time_range)
            disp([num2str(tt) '/' num2str(length(time_range))])
            itime = time_range(tt);
            
            data_type_A = 4; % Get Baltic+ NOM (v1.0)
            [A] = rd_smos_BEC_Baltic_FUN(itime,ibasin,data_type_A);
            data_type_B = 5; % Get Baltic+ NOM (v1.0)
            [B] = rd_smos_BEC_Baltic_FUN(itime,ibasin,data_type_B);
            
            
            sss_nm = A.sss;
            ssse_nm = A.sss_error;
            sss_nm_irange = A.sss_irange;
            sss_error_nm_irange = A.sss_error_irange;
            
            sss_ns = B.sss;
            ssse_ns = B.sss_error;
            sss_ns_irange = B.sss_irange;
            sss_error_ns_irange = B.sss_error_irange;
            
            lon_smos = A.lon_irange;
            lat_smos = A.lat_irange;
            
            % Set dimension loop-output (SSS and SSS_error)
            if tt == 1
                [a,b,c] = size(sss_nm_irange);
                z = length(time_range);
                
                % SSS-SMOS at whole Baltic
                SSS_NM = nan(a,b,z);
                SSSe_NM = nan(a,b,z);
                
                SSS_NS = nan(a,b,z);
                SSSe_NS = nan(a,b,z);
                
                
                % SSS-SMOS in each Baltic-region
                SSS_NM_irange = nan(a,b,c,z);
                SSSe_NM_irange = nan(a,b,c,z);
                
                SSS_NS_irange = nan(a,b,c);
                SSSe_NS_irange = nan(a,b,c);
            end
            
            SSS_NM (:,:,tt) = sss_nm;
            SSSe_NM (:,:,tt) = ssse_nm;
            
            SSS_NS (:,:,tt) = sss_ns;
            SSSe_NS (:,:,tt) = ssse_ns;
            
            
            SSS_NM_irange(:,:,:,tt) = sss_nm_irange;
            SSSe_NM_irange(:,:,:,tt) = sss_error_nm_irange;
            
            SSS_NS_irange(:,:,:,tt) = sss_ns_irange;
            SSSe_NS_irange(:,:,:,tt) = sss_error_ns_irange;
        end
        
        SSS_NM = permute(SSS_NM_irange,[1,2,4,3]);
        SSSe_NM = permute(SSSe_NM_irange,[1,2,4,3]);
        SSS_NS = permute(SSS_NS_irange,[1,2,4,3]);
        SSSe_NS = permute(SSSe_NS_irange,[1,2,4,3]);
        
        SSS_NM_irange = permute(SSS_NM_irange,[1,2,4,3]);
        SSSe_NM_irange = permute(SSSe_NM_irange,[1,2,4,3]);
        SSS_NS_irange = permute(SSS_NS_irange,[1,2,4,3]);
        SSSe_NS_irange = permute(SSSe_NS_irange,[1,2,4,3]);
        
        save(fn_out,vars_save{:})
        
    elseif fn_exist == 2
        clc
        load (fn_out)
        disp('SSS-SMOS in each region data: ')
        
        whos (vars_save{:})
    end
    
    ind1 = SSS_NM_irange<0;
    ind2 = SSS_NS_irange<0;
    
    SSS_NM_irange(ind1) = [];
    SSSe_NM_irange(ind1) = [];
    
    SSS_NS_irange(ind2) = [];
    SSSe_NS_irange(ind2) = [];
    
    %% [4] SSS-SMOS STATS
    
    %% [4.1] SSS-SMOS Monthly stats
    iMONTHS = 1:12;
    iYEARS = iyear1:iyear2;
    
    nYEARS = length(iYEARS);
    nMONTHS = length(iMONTHS);
    
    
    % Month-numbers in itime (SMOS)
    itime_vec = datevec (time_range);
    
    [a,b,c,d] = size(SSS_NM_irange);
    
    c = nYEARS * nMONTHS;
    
    % Whole Baltic SSS-SMOS
    SSS_NM_mean_month = nan(a,b,c);
    SSSe_NM_mean_month = nan(a,b,c);
    SSS_NS_mean_month = nan(a,b,c);
    SSSe_NS_mean_month = nan(a,b,c);
    
    
    % Balic-Regional SSS-SMOS
    AA = nan(a,b,c,d);
    
    SSS_NM_irange_mean_month = AA;
    SSSe_NM_irange_mean_month = AA;
    
    SSS_NS_irange_mean_month = AA;
    SSSe_NS_irange_mean_month = AA;
    
    ind_position = 1:c;
    
    iMONTHS2 = repmat(iMONTHS,1,nYEARS);
    iYEARS2 = repmat(iYEARS,nMONTHS,1);
    
    iYEARS = vertcat(iYEARS2(:,1),iYEARS2(:,2),iYEARS2(:,3));
    iMONTHS = iMONTHS2;
    
    clear iYEARS2 iMONTHS2
    
    % =================
    % Loop Compute: SSS-SMOS month averages at basin and at subregion scale
    % e.g [Basin]: SSS_NM and [sub-basin]: SSS_NM_irange
    %
    for nn = 1: c
        this_year = iYEARS(nn);
        this_month = iMONTHS(nn);
        
        % index month-number and make monthly average
        ind  = find(itime_vec(:,1) == this_year & itime_vec (:,2) == this_month);
        
        if ~isempty(ind)
            SSS_NM_mean_month_alfa = nanmean(SSS_NM(:,:,ind),3);
            SSS_NM_mean_month (:,:,nn) = SSS_NM_mean_month_alfa;
            
            SSSe_NM_mean_month_alfa = nanmean(SSSe_NM(:,:,ind),3);
            SSSe_NM_mean_month (:,:,nn) = SSSe_NM_mean_month_alfa;
            
            SSS_NS_mean_month_alfa = nanmean(SSS_NS(:,:,ind),3);
            SSS_NS_mean_month (:,:,nn) = SSS_NS_mean_month_alfa;
            
            SSSe_NS_mean_month_alfa = nanmean(SSSe_NS(:,:,ind),3);
            SSSe_NS_mean_month (:,:,nn) = SSSe_NS_mean_month_alfa;
            
            
            
            SSS_NM_irange_mean_month_alfa = nanmean(SSS_NM_irange(:,:,ind,:),3);
            SSS_NM_irange_mean_month (:,:,nn,:) = SSS_NM_irange_mean_month_alfa;
            
            SSSe_NM_irange_mean_month_alfa = nanmean(SSSe_NM_irange(:,:,ind,:),3);
            SSSe_NM_irange_mean_month (:,:,nn,:) = SSSe_NM_irange_mean_month_alfa;
            
            SSS_NS_irange_mean_month_alfa = nanmean(SSS_NS_irange(:,:,ind,:),3);
            SSS_NS_irange_mean_month (:,:,nn,:) = SSS_NS_irange_mean_month_alfa;
            
            SSSe_NS_irange_mean_month_alfa = nanmean(SSS_NS_irange(:,:,ind,:),3);
            SSSe_NS_irange_mean_month (:,:,nn,:) = SSS_NS_irange_mean_month_alfa;
            
            clear *_alfa
            
        end
    end
    %
    % =================
    
    
    %% [4.2] SSS-SMOS Regional stats
    [a,b,c,d] = size(SSS_NM_irange);
    
    SSS_NM_irange2 = reshape(SSS_NM_irange,a*b,c,d);
    SSSe_NM_irange2 = reshape(SSSe_NM_irange,a*b,c,d);
    
    SSS_NS_irange2 = reshape(SSS_NS_irange,a*b,c,d);
    SSSe_NS_irange2 = reshape(SSSe_NS_irange,a*b,c,d);
    
    SSS_NM_irange_mean = squeeze(nanmean(SSS_NM_irange2,1));
    SSSe_NM_irange_mean = squeeze(nanmean(SSSe_NM_irange2,1));
    
    SSS_NS_irange_mean = squeeze(nanmean(SSS_NS_irange2,1));
    SSSe_NS_irange_mean = squeeze(nanmean(SSSe_NS_irange2,1));
    
    SSS_NM_irange_std = squeeze(nanstd(SSS_NM_irange2,1));
    SSSe_NM_irange_std = squeeze(nanstd(SSSe_NM_irange2,1));
    
    SSS_NS_irange_std = squeeze(nanstd(SSS_NS_irange2,1));
    SSSe_NS_irange_std = squeeze(nanstd(SSSe_NS_irange2,1));
    
    
    %% [5] SSS-SMOS Time Series
    
    %% [5.1] SSS Time series Raw-series
    
    for nn = 1:length(SREG_str)
        
        X1 = time_range';
        Y1 = SSS_NM_irange_mean(:,nn);
        Z1 = SSS_NM_irange_std(:,nn);
        
        X2 = time_range';
        Y2 = SSS_NS_irange_mean(:,nn);
        Z2 = SSS_NS_irange_std(:,nn);
        
        
        ind1 = isnan(Y1) | Y1<0 | Z1<0 |isnan(Z1);
        ind2 = isnan(Y2) | Y2<0 | Z2<0 |isnan(Z2);
        
        X1(ind1) = [];
        Y1(ind1) = [];
        Z1(ind1) = [];
        
        X2(ind2) = [];
        Y2(ind2) = [];
        Z2(ind2) = [];
        
        % ====================================
        figure(nn); clf
        %set(gcf,'DefaultAxesFontSize',24);
        
        
        h1 = boundedline(X1,Y1,Z1,'r-',...
            'alpha','orientation','vert','transparency',0.3); hold on
        set(h1,'linewidth',2);
        
        h2 = boundedline(X2,Y2,Z2,'y-',...
            'alpha','orientation','vert','transparency',0.3); hold on
        set(h2,'linewidth',2);
        
        h3 = boundedline(TIME_model,SSS_model_mean_sreg(:,nn),SSS_model_std_sreg(:,nn), 'b-',...
            'alpha','orientation','vert','transparency',0.3); hold on
        set(h3,'linewidth',2);
        
        lg = legend([h1(1);h2(1);h3(1)],...
            {' SSS-BEC_{NM}';' SSS-BEC_{NS}' ; 'SSS_{NEMO}'},'fontsize',14);
        set(lg,'Location','northwest')
        
        box on
        set(gca,'XMinorTick','on','YMinorTick','on')
        
        grid minor
        axis tight
        
        
        set(gca,'xtick',TIME_model(1:2:end))
        datetick('x','mm/yy','keepticks')
        title(['Regional averaged (mean) Baltic ' SREG_str2{nn}],'fontsize',24);
        xlabel('time (months)','fontsize',24);
        ylabel([ 'SSS (psu)' ], 'fontsize',24);
        
        fg_name = [basin_str '_MEAN_NEMOvsSAT_' SREG_str2{nn} '_'...
            time_str1 '_' time_str2];
        
        % Save figure  - output -
        folder_this = [folder_figs 'NEMOvsSAT/TIMESERIES/MEAN/'];
        
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
    
    %% [5.2] Basin wide SSS Time series differences (dSSS = Satellite minus Model)
    
    dSSS_NM = SSS_NM_mean_month - SSS_model;
    dSSS_NS = SSS_NS_mean_month - SSS_model;
    
    [a,b,c] = size(dSSS_NM);
    dSSS_NM_alfa = reshape(dSSS_NM,a*b,c);
    dSSS_NS_alfa = reshape(dSSS_NS,a*b,c);
    
    dSSS_NM_mean = nanmean(dSSS_NM_alfa,1);
    dSSS_NM_std = nanstd(dSSS_NM_alfa,1);
    
    dSSS_NS_mean = nanmean(dSSS_NS_alfa,1);
    dSSS_NS_std = nanstd(dSSS_NS_alfa,1);
    
    clear *_alfa
    
    X1 = TIME_model';
    Y1 = dSSS_NM_mean;
    Z1 = dSSS_NS_std;
    
    X2 = TIME_model';
    Y2 = dSSS_NS_mean;
    Z2 = dSSS_NS_std;
    
    
    ind1 = isnan(Y1) | isnan(Z1);
    ind2 = isnan(Y2) | isnan(Z2);
    
    X1(ind1) = [];
    Y1(ind1) = [];
    Z1(ind1) = [];
    
    X2(ind2) = [];
    Y2(ind2) = [];
    Z2(ind2) = [];
    
    
    % ====================================
    figure; clf
    %set(gcf,'DefaultAxesFontSize',24);
    
    h1 = boundedline(X1,Y1,Z1,'r-',...
        'alpha','orientation','vert','transparency',0.3); hold on
    set(h1,'linewidth',2);
    
    h2 = boundedline(X2,Y2,Z2,'y-',...
        'alpha','orientation','vert','transparency',0.3); hold on
    set(h2,'linewidth',2);
    
    
    lg = legend([h1(1);h2(1)],...
        {' SSS-BEC_{NM} - SSS_{NEMO}';' SSS-BEC_{NS} - SSS_{NEMO}' },'fontsize',14);
    set(lg,'Location','northwest')
    
    box on
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    grid minor
    axis tight
    
    
    set(gca,'xtick',TIME_model(1:2:end))
    datetick('x','mm/yy','keepticks')
    title(['Regional averaged (mean) Baltic [ALL]'],'fontsize',24);
    xlabel('time (months)','fontsize',24);
    ylabel([ 'SSS (psu)' ], 'fontsize',24);
    
    fg_name = [basin_str '_dSSS_NEMOvsSAT_ALL_'...
        time_str1 '_' time_str2];
    
    % Save figure  - output -
    folder_this = [folder_figs 'NEMOvsSAT/TIMESERIES/dSSS/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name '.' fg_format];
    
    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
    
    
    %% [5.3] SUB-Basin SSS Time series differences (dSSS = Satellite minus Model)
    
    [a,b,c,d] = size(SSS_NM_irange_mean_month);
    
    SSS_NM_irange_mean_month2 = reshape(SSS_NM_irange_mean_month,a*b,c,d);
    SSSe_NM_irange_mean_month2 = reshape(SSSe_NM_irange_mean_month,a*b,c,d);
    
    SSS_NS_irange_mean_month2 = reshape(SSS_NS_irange_mean_month,a*b,c,d);
    SSSe_NS_irange_mean_month2 = reshape(SSSe_NS_irange_mean_month,a*b,c,d);
    
    
    dSSS_NM_irange = SSS_NM_irange_mean_month2 - SSS_model_irange;
    dSSS_NS_irange = SSS_NS_irange_mean_month2 - SSS_model_irange;
    
    
    dSSS_NM_irange_mean = squeeze(nanmean(dSSS_NM_irange,1));
    dSSS_NM_irange_std = squeeze(nanstd(dSSS_NM_irange,1));
    
    dSSS_NS_irange_mean = squeeze(nanmean(dSSS_NS_irange,1));
    dSSS_NS_irange_std = squeeze(nanstd(dSSS_NS_irange,1));
    
    %% Plot dSSS (Sat - Model)
    for nn = 1:length(SREG_str)
        
        
        X1 = TIME_model';
        Y1 = dSSS_NM_irange_mean(:,nn);
        Z1 = dSSS_NS_irange_std(:,nn);
        
        X2 = TIME_model';
        Y2 = dSSS_NS_irange_mean(:,nn);
        Z2 = dSSS_NS_irange_std(:,nn);
        
        
        ind1 = isnan(Y1) | isnan(Z1);
        ind2 = isnan(Y2) | isnan(Z2);
        
        X1(ind1) = [];
        Y1(ind1) = [];
        Z1(ind1) = [];
        
        X2(ind2) = [];
        Y2(ind2) = [];
        Z2(ind2) = [];
        
        % ====================================
        figure; clf
        %set(gcf,'DefaultAxesFontSize',24);
        
        h1 = boundedline(X1,Y1,Z1,'r-',...
            'alpha','orientation','vert','transparency',0.3); hold on
        set(h1,'linewidth',2);
        
        h2 = boundedline(X2,Y2,Z2,'y-',...
            'alpha','orientation','vert','transparency',0.3); hold on
        set(h2,'linewidth',2);
        
        plot(X2,Y1*.0,'k-','linewidth',2)
        
        axis tight
        
        lg = legend([h1(1);h2(1)],...
            {' SSS-BEC_{NM} - SSS_{NEMO}';' SSS-BEC_{NS} - SSS_{NEMO}' },'fontsize',14);
        set(lg,'Location','northwest')
        
        box on
        set(gca,'XMinorTick','on','YMinorTick','on')
        
        grid minor
        
        
        set(gca,'xtick',TIME_model(1:2:end))
        datetick('x','mm/yy','keepticks')
        title(['Regional averaged (mean) dSSS Baltic ' SREG_str2{nn}],'fontsize',24);
        xlabel('time (months)','fontsize',24);
        ylabel([ 'dSSS (psu)' ], 'fontsize',24);
        
        fg_name = [basin_str '_dSSS_NEMOvsSAT_' SREG_str2{nn} '_'...
            time_str1 '_' time_str2];
        
        % Save figure  - output -
        folder_this = [folder_figs 'NEMOvsSAT/TIMESERIES/dSSS/'];
        
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
    
    
    % BOXPLOTS
    %% [1] BOXPLOT ALL-BALTIC - BOXplot need to be updated <if needed>
    % iparam_str_all = {'SSS','SST','MLD'};
    
    % figure; clf; hold on;
    % set(gcf,'DefaultAxesFontSize',24);
    % set(gca,'TickLabelInterpreter','tex');
    %
    % bplot = boxplot...
    %     ([SSS_model(:) SSS_NM_irange2]);
    %
    % ylim([0 10]);
    %
    % ylabel('S')
    %
    % h = findobj(gca, 'type', 'text');
    % set(h,'FontSize',25);
    % set(h,'VerticalAlignment', 'middle');
    % set(gca,'FontSize',25);
    %
    % title({['Satellite SSS around Argo float (plat ID:'...
    %     num2str(platform(nn))...
    %     ' , PROF#' num2str(ID(nn)) ')'];...
    %     ['Distance: ' num2str(r) ' km']});
    %
    % fg_name = [fn '_SAT-SSS_BOXPLOT_PLAT'...
    %     num2str(platform(nn))...
    %     '_PROF#' num2str(ID(nn))...
    %     '_R' num2str(r) 'KM.'...
    %     fg_format];
    %
    % % Save figure  - output -
    % folder_this = [folder_figs 'BOXPLOTS/'];
    %
    % if fg_save == 1
    %     foldercheck_raf(folder_this); %! make folder_figs
    % end
    %
    % fg_name = [folder_this fg_name];
    %
    %
    % fg_exist = exist(fg_name,'file'); % check fn existence
    %
    % if fg_save == 1 && fg_exist == 0
    %     save_raf(gcf,fg_name,fg_format); close
    % end
    %
    
    
    
end
