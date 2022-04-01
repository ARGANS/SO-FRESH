% Syntax: Baltic_argo_BEC_analyses.m (script)
%
% Description
% Analyses of Argo floats as provided by BEC.
% Input: matlab files created from the original netCDF files (use
% ar_rd_BEC.m)
%
% This script compares Satellite data againts the Argo (in situ) data.
% Furthermore, this script does compare Argo againts regional model data (NEMO).
% The model data is used in the processing of the new Baltic+ SSS data set.
%
% This script does the following:
% [1] Colocations Argo-2-satellite and Argo-2-Model grid/time [Baltic_argo_colocations.m]
% [2] Comparisons Argo Vs Satellite and Argo-2-Model
% [3] Baltic Study [2011-2013]:
% *[3.1] BOXPLOTS
% ** [3.1.1] BOXPPLOT: SSS at float location (i.e. SSS-satellite and SSS-model)
% *[3.2] PROFILES: Argo Profiles (T and S) -- study column structure
% *[3.3] TIMESERIES
% ** [3.3.1] SSS Timeseries at float location
% ** [3.3.2] dSSS Timeseries at float location (dSSS = Argo - QD, QD: Query dataset)
%
%
% Data source:
% BEC products for validation
%
%
% current version: v1r3 (2020/02/11)
%
% History
%
% [3] Change order of script by:
%     --> 1/ Make daily plots Argo and SSS-satellite;
%     --> 2/ Save daily plots (to study TOP surface melting ice Vs. SSS)
%
% [2] Creation of snipet (functions) to make Argo vertical interpolations
% [1] Creation of this script [20191114]
%
%
% =========================================================================
% Author: rcatany
%
% =========================================================================

%clc;
clear;
close all

% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance platform (in Km)
depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m


fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/argo/argo_mat/']);

% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'ARGO_MISSING_20200121.txt'];

folder_figs = '/Volumes/Rogue/scratch/Validation/';

% ==========================
iyear = 2013;
imonth = 1:12;


% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file

ibasin = 9; % Basin number: 7 (Arctic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

folder_figs = [folder_figs basin_str '/Argo/Argo-BEC/'];

if fg_save == 1
    foldercheck_raf(folder_figs); %! make folder_figs
end


%% [1] Colocation Argo-2-Satellite
run Baltic_argo_BEC_colocations.m % (snippet)


%% [1.1] Compute Argo at 10 m and the Argo average surface-2-10m
nprof = length(platform_irange);

[zlevels,~] = size(PRES_intp);

SALT10_mean = nan(1,nprof);
SALT10_std = nan(1,nprof);
SALT10_point = nan(1,nprof);

for nn = 1:nprof
    %index the average top 10m and find the value at 10m
    ind_mean = find(PRES_intp(:,nn)<=30,zlevels,'first');
    ind_point = find(PRES_intp(:,nn)==10,zlevels,'first');
    
    SALT_beta = SALT_intp(ind_mean,nn);
    SALT10_mean(1,nn) = nanmean(SALT_beta);
    SALT10_std (1,nn) = nanstd(SALT_beta);
    
    SALT10_point_beta = SALT_intp(ind_point,nn);
    
    while isnan(SALT10_point_beta)
        ind_point = ind_point+1;
        SALT10_point_beta = SALT_intp(ind_point,nn);
    end
    
    SALT10_point(1,nn) = SALT10_point_beta;
    
end; clear nn *beta*


%% [2] PLOTS: Argo Vs Satellite
%% [2.1] Plot: Argo geolocation
% Plots in one map locations of same platform

make_plot = 1; % flag to make plot [1], or not [0]
if make_plot == 1
    
    NPLATS = unique(platform_irange);
    
    for nn = 1: length(NPLATS)
        
        % snippet to plot map and zoom-in to Argo location
        plot_map_argo_BEC
        
    end
end



%% [2.2] PLOT Argo profiles
%% [2.2.1] plot raw Vs interp T and S of a given Profile (nprof)

nprof = length(platform_irange);

make_plots = 0;

if make_plots == 1
    for nn = 1:nprof
        
        TEMP_beta1 = TEMP_irange(:,nn);
        SALT_beta1 = SALT_irange(:,nn);
        PRES_beta1 = PRES_irange(:,nn);
        
        TEMP_beta2 = TEMP_intp(:,nn);
        SALT_beta2 = SALT_intp(:,nn);
        PRES_beta2 = PRES_intp(:,nn);
        
        iTIME = [time_range(nn)];
        
        str1 = 'RAW';
        str2 = 'intp';
        
        
        make_plot = 1; % flag to make plot [1], or not [0]
        if make_plot == 1
            
            figure(4); clf; hold on;
            set(gcf,'DefaultAxesFontSize',18);
            
            plot_profile_argo_BEC(TEMP_beta1,SALT_beta1,PRES_beta1,...
                TEMP_beta2,SALT_beta2,PRES_beta2,platform_irange(nn),...
                iTIME,str1,str2)
            
            % Save figure  - output -
            fg_name = ['argo_' datestr(iTIME,'yyyymmdd') '_RAWvsINTP_PLAT'...
                num2str(platform_irange(nn))...
                '_PROF#' num2str(ID_irange(nn)) ...
                '.'  fg_format];
            
            % Save figure  - output -
            folder_this = [folder_figs 'RAWvsINT/'];
            
            if fg_save == 1
                foldercheck_raf(folder_this); %! make folder_figs
            end
            
            fg_name = [folder_this fg_name];
            
            % check fn existence
            fg_exist = exist(fg_name,'file');
            if fg_save == 1 && fg_exist == 0
                save_raf(gcf,fg_name,fg_format); close
            end
        end
        
    end; clear nn *_beta*
    
end



%% [2.2.2] COLORPLOT all-Profiles by the same Argo platform
NPLATS = unique(platform_irange);

for nn = 1: length(NPLATS)
    
    this_platform = NPLATS(nn);
    
    ind = find(platform_irange == this_platform);
    
    time_ar = time_range(ind);
    
    time_str1 = datestr(time_ar(1),'yyyymmdd');
    time_str2 = datestr(time_ar(end),'yyyymmdd');
    
    TEMP_beta1 = TEMP_irange(:,ind);
    SALT_beta1 = SALT_irange(:,ind);
    PRES_beta1 = PRES_irange(:,ind);
    
    TEMP_beta2 = TEMP_intp(:,ind);
    SALT_beta2 = SALT_intp(:,ind);
    PRES_beta2 = PRES_intp(:,ind);
    
    colorplot_argo_platform_profiles(...
        this_platform,time_ar,...
        TEMP_beta1,SALT_beta1,PRES_beta1,...
        TEMP_beta2,SALT_beta2,PRES_beta2);
    
    % Save figure  - output -
    fg_name = ['argo_' num2str(this_platform) '_COLORPLOT_'...
        time_str1 '_' time_str2 '.'  fg_format];
    
    % Save figure  - output -
    folder_this = [folder_figs 'COLORPLOT/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name];
    
    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
end; clear ind nn *beta*


%% [2.2.3] More than one profile-platform, make mean and compare all profiles
nprof = length(platform_irange);

NPLATS = unique(platform_irange);

make_plot = 1;

if make_plot == 1
    
    for nn = 1:length(NPLATS)
        
        this_platform = NPLATS(nn);
        
        ind = find(platform_irange == this_platform);
        
        time_ar = time_range(ind);
        
        time_str1 = datestr(time_ar(1),'yyyymmdd');
        time_str2 = datestr(time_ar(end),'yyyymmdd');
        
        
        TEMP_beta1 = TEMP_intp(:,ind);
        SALT_beta1 = SALT_intp(:,ind);
        PRES_beta1 = PRES_intp(:,ind);
        
        TEMP_beta2 = TEMP_intp(:,ind);
        SALT_beta2 = SALT_intp(:,ind);
        PRES_beta2 = PRES_intp(:,ind);
        
        
        str1 = 'all-profiles';
        str2 = 'mean';
        
        % ===========
        
        boundedplot_flag = 1;
        
        figure(4); clf; hold on
        set(gcf,'DefaultAxesFontSize',18);
        
        plot_profile_argo_BEC(TEMP_beta1,SALT_beta1,PRES_beta1,...
            TEMP_beta2,SALT_beta2,PRES_beta2,this_platform,...
            time_ar,str1,str2,boundedplot_flag)
        
        % Save figure  - output -
        fg_name = ['ar_' num2str(this_platform) '_'...
            time_str1 '_' time_str2 '_BOUNDEDPLOT_PLAT'...
            '.'  fg_format];
        
        % Save figure  - output -
        folder_this = [folder_figs 'BOUNDEDPLOT/'];
        
        if fg_save == 1
            foldercheck_raf(folder_this); %! make folder_figs
        end
        
        fg_name = [folder_this fg_name];
        
        % check fn existence
        fg_exist = exist(fg_name,'file');
        
        if fg_save == 1 && fg_exist == 0
            save_raf(gcf,fg_name,fg_format); close
        end
        
    end
    
end; clear nn *_beta*



%% [2.4] BOXPLOTS
%% [2.4.1] boxplots: Argo data versus Satellite data (at float location)
% shall we use the median of SSS to be compared againts

make_plot = 1; % flag to make plot [1], or not [0]
if make_plot == 1
    % string product name
    PRODS = {...
        'SSS_{NS}','SSS_{NM}','SSS_{gl}',...
        'SSS_{NEMO}'};
    
    
    NPLATS = unique(platform_irange);
    
    for nn = 1:length(NPLATS)
        
        this_platform = NPLATS(nn);
        
        ind = find(platform_irange == this_platform);
        
        nprof = length(ind);
        
        time_ar = time_range(ind);
        
        time_str1 = datestr(time_ar(1),'yyyymmdd');
        time_str2 = datestr(time_ar(end),'yyyymmdd');
        
        
        X_beta = sss_ns_irange(:,ind);
        Y_beta = sss_nm_irange(:,ind);
        Z_beta = sss_gl_irange(:,ind);
        T_beta = sss_model_irange(:,ind);
        
        X_beta = X_beta(:);
        Y_beta = Y_beta(:);
        Z_beta = Z_beta(:);
        T_beta = T_beta(:);
        
        X_beta (X_beta == 0) = NaN;
        Y_beta (Y_beta == 0) = NaN;
        Z_beta (Z_beta == 0) = NaN;
        T_beta (T_beta == 0) = NaN;
        
        figure; clf; hold on;
        set(gcf,'DefaultAxesFontSize',24);
        
        bplot = boxplot([X_beta,Y_beta,Z_beta,T_beta],'labels',PRODS);
        
        bp = gca;
        bp.XAxis.TickLabelInterpreter = 'tex';
        
        ylim([0 15]);
        
        ylabel('SSS')
        
        grid on
        grid minor
        
        h = findobj(gca, 'type', 'text');
        set(h,'FontSize',25);
        set(h,'VerticalAlignment', 'middle');
        set(gca,'FontSize',25);
        
        title({...
            [basin_str ' Satellite SSS around Argo float '];...
            ['ID: ' num2str(this_platform) ...
            ' (' num2str(nprof) ' profiles)'...
            ', Distance: ' num2str(r) ' km'];...
            [time_str1 '-' time_str2]});
        
        fg_name = ['ar_PLAT#' num2str(platform_irange(nn))...
            '_SAT-SSS_BOXPLOT_'...
            time_str1 '_' time_str2...
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
        
    end; clear ind *beta* *exist*
    
end

%% [2.4.2] boxplots: ALL Argo data versus ALL Satellite data

make_plot = 1; % flag to make plot [1], or not [0]
if make_plot == 1
    % string product name
    PRODS = {...
        'SSS_{NS}','SSS_{NM}','SSS_{gl}',...
        'SSS_{NEMO}'};
    
    nprof = length(platform_irange);
    
    time_ar = time_range(:);
    
    time_str1 = datestr(time_ar(1),'yyyymmdd');
    time_str2 = datestr(time_ar(end),'yyyymmdd');
    
    
    X_beta = sss_ns_irange(:);
    Y_beta = sss_nm_irange(:);
    Z_beta = sss_gl_irange(:);
    T_beta = sss_model_irange(:);
    
    X_beta = X_beta(:);
    Y_beta = Y_beta(:);
    Z_beta = Z_beta(:);
    T_beta = T_beta(:);
    
    X_beta (X_beta == 0) = NaN;
    Y_beta (Y_beta == 0) = NaN;
    Z_beta (Z_beta == 0) = NaN;
    T_beta (T_beta == 0) = NaN;
    
    figure; clf; hold on;
    set(gcf,'DefaultAxesFontSize',24);
    
    bplot = boxplot([X_beta,Y_beta,Z_beta,T_beta],'labels',PRODS);
    
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'tex';
    
    ylim([0 15]);
    
    ylabel('SSS')
    
    grid on
    grid minor
    
    h = findobj(gca, 'type', 'text');
    set(h,'FontSize',25);
    set(h,'VerticalAlignment', 'middle');
    set(gca,'FontSize',25);
    
    title({...
        [basin_str ' Satellite SSS around Argo float '];...
        ['ALL-platforms  (' num2str(nprof) ' profiles)'...
        ', Distance: ' num2str(r) ' km'];...
        [time_str1 '-' time_str2]});
    
    fg_name = ['ar_PLAT#ALL'...
        '_SAT-SSS_BOXPLOT_'...
        time_str1 '_' time_str2...
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
    
end; clear ind *beta* *exist*



%% [2.5] TIMESERIES
%% [2.5.1] TIMESERIES: SSS at float location
% Compare each product to Salinity Argo 10-m and to Argo [0-10m] average

% string product name
PRODS = {...
    ' SSS_{NS}',' SSS_{NM}',' SSS_{gl}',...
    ' SSS_{NEMO}', ' S_{Argo-10}',' S_{Argo-10*}'};

make_plot = 1;

if make_plot == 1
    
    % Get more colors: type 'uisetcolor'
    color_lines = [...
        0    0.4471    0.7412;...
        1.0000    1.0000    0.0667;...
        0.3922    0.8314    0.0745;...
        0.7176    0.2745    1.0000;...
        0.9294    0.6941    0.1255;...
        1.0000    0.0745    0.6510];
    
    
    NPLATS = unique(platform_irange);
    
    for nn = 1: length(NPLATS)
        
        this_platform = NPLATS(nn);
        
        ind = find(platform_irange == this_platform);
        
        nprof = length(ind);
        
        X = time_range;
        time_line = X(1):20:X(end);
        time_line_diff = diff(time_line);
        
        itime_start_str = datestr(X(1),'yyyymmdd');
        itime_end_str = datestr(X(end),'yyyymmdd');
        
        Y_beta1 = sss_ns_irange_mn(ind);
        Y_beta2 = sss_nm_irange_mn(ind);
        Y_beta3 = sss_gl_irange_mn(ind);
        Y_beta4 = sss_model_irange_mn(ind);
        Y_beta5 = SALT10_point(ind);
        Y_beta6 = SALT10_mean(ind);
        
        Y_beta1 (Y_beta1 == 0) = NaN;
        Y_beta2 (Y_beta2 == 0) = NaN;
        Y_beta3 (Y_beta3 == 0) = NaN;
        Y_beta4 (Y_beta4 == 0) = NaN;
        Y_beta5 (Y_beta5 == 0) = NaN;
        Y_beta6 (Y_beta6 == 0) = NaN;
        
        % ==========================
        figure
        [hl, hp] = boundedline(X,SALT10_mean(ind),SALT10_std(ind),'-k',...
            'transparency',0.5);
        
        hold on
        
        h1 = plot(X,Y_beta1,'-o','linewidth',2,'color',color_lines(1,:),...
            'MarkerFaceColor',color_lines(1,:));
        h2 = plot(X,Y_beta2,'-o','linewidth',2,'color',color_lines(2,:),...
            'MarkerFaceColor',color_lines(2,:));
        h3 = plot(X,Y_beta3,'-o','linewidth',2,'color',color_lines(3,:),...
            'MarkerFaceColor',color_lines(3,:));
        h4 = plot(X,Y_beta4,'-o','linewidth',2,'color',color_lines(4,:),...
            'MarkerFaceColor',color_lines(4,:));
        h5 = plot(X,Y_beta5,'-k','linewidth',2);
        h6 = plot(X,Y_beta6,'-r','linewidth',2);
        
        xlabel ('Time (days)','fontsize',20)
        ylabel ('SSS (psu)','fontsize',20)
        
        grid on
        grid minor
        box on
        
        h = findobj(gca, 'type', 'text');
        set(h,'FontSize',25);
        set(h,'VerticalAlignment', 'middle');
        set(gca,'FontSize',25);
        
        xticks(time_line)
        xticklabels(datestr(time_range,'ddmmyy'))
        xtickangle(45)
        
        lg = legend([h1,h2,h3,h4,h5,h6],PRODS,...
            'fontsize',18,'location','NorthWest');
        
        title({...
            [basin_str ' Satellite SSS around Argo float '];...
            ['ID: ' num2str(this_platform) ' (' num2str(nprof) ' profiles)'];...
            [itime_start_str '-' itime_end_str]
            ['Distance: ' num2str(r) ' km']});
        
        fg_name = [...
            'ar_PLAT#' num2str(this_platform)...
            '_ARG2SAT_SSS_TIMESERIES_'...
            itime_start_str '_' itime_end_str...
            '_R' num2str(r) 'KM.'...
            fg_format];
        
        % Save figure  - output -
        folder_this = [folder_figs 'TIMESERIES/SSS/'];
        
        if fg_save == 1
            foldercheck_raf(folder_this); %! make folder_figs
        end
        
        fg_name = [folder_this fg_name];
        
        % check fn existence
        fg_exist = exist(fg_name,'file');
        
        if fg_save == 1 && fg_exist == 0
            save_raf(gcf,fg_name,fg_format); close
        end
        
    end
    
end
%% [2.5.2] TIMESERIES: ALL-SSS at ALL-float location
% Compare each product to Salinity Argo 10-m and to Argo [0-10m] average

% string product name
PRODS = {...
    ' SSS_{NS}',' SSS_{NM}',' SSS_{gl}',...
    ' SSS_{NEMO}', ' S_{Argo-10}',' S_{Argo-10*}'};

make_plot = 1;

if make_plot == 1
    % Get more colors: type 'uisetcolor'
    color_lines = [...
        0    0.4471    0.7412;...
        1.0000    1.0000    0.0667;...
        0.3922    0.8314    0.0745;...
        0.7176    0.2745    1.0000;...
        0.9294    0.6941    0.1255;...
        1.0000    0.0745    0.6510];
    
    nprof = length(platform_irange);
    
    X = time_range;
    time_line = X(1):20:X(end);
    time_line_diff = diff(time_line);
    
    time_str1 = datestr(X(1),'yyyymmdd');
    time_str2 = datestr(X(end),'yyyymmdd');
    
    Y_beta1 = sss_ns_irange_mn(:);
    Y_beta2 = sss_nm_irange_mn(:);
    Y_beta3 = sss_gl_irange_mn(:);
    Y_beta4 = sss_model_irange_mn(:);
    Y_beta5 = SALT10_point(:);
    Y_beta6 = SALT10_mean(:);
    
    Y_beta1 (Y_beta1 == 0) = NaN;
    Y_beta2 (Y_beta2 == 0) = NaN;
    Y_beta3 (Y_beta3 == 0) = NaN;
    Y_beta4 (Y_beta4 == 0) = NaN;
    Y_beta5 (Y_beta5 == 0) = NaN;
    Y_beta6 (Y_beta6 == 0) = NaN;
    
    % ==========================
    figure
    [hl, hp] = boundedline(X,SALT10_mean(:),SALT10_std(:),'-k',...
        'transparency',0.5);
    
    hold on
    
    h1 = plot(X,Y_beta1,'-o','linewidth',2,'color',color_lines(1,:),...
        'MarkerFaceColor',color_lines(1,:));
    h2 = plot(X,Y_beta2,'-o','linewidth',2,'color',color_lines(2,:),...
        'MarkerFaceColor',color_lines(2,:));
    h3 = plot(X,Y_beta3,'-o','linewidth',2,'color',color_lines(3,:),...
        'MarkerFaceColor',color_lines(3,:));
    h4 = plot(X,Y_beta4,'-o','linewidth',2,'color',color_lines(4,:),...
        'MarkerFaceColor',color_lines(4,:));
    h5 = plot(X,Y_beta5,'-k','linewidth',2);
    h6 = plot(X,Y_beta6,'-r','linewidth',2);
    
    xlabel ('Time (days)','fontsize',20)
    ylabel ('SSS (psu)','fontsize',20)
    
    grid on
    grid minor
    box on
    
    h = findobj(gca, 'type', 'text');
    set(h,'FontSize',25);
    set(h,'VerticalAlignment', 'middle');
    set(gca,'FontSize',25);
    
    xticks(time_line)
    xticklabels(datestr(time_range,'ddmmyy'))
    xtickangle(45)
    
    lg = legend([h1,h2,h3,h4,h5,h6],PRODS,...
        'fontsize',18,'location','NorthWest');
    
    title({...
        [basin_str ' Satellite SSS around Argo float '];...
        ['ALL-platforms (' num2str(nprof) ' profiles)'];...
        [time_str1 '-' time_str2]
        ['Distance: ' num2str(r) ' km']});
    
    fg_name = [...
        'ar_PLAT#ALL' ...
        '_ARG2SAT_SSS_TIMESERIES_'...
        time_str1 '_' time_str2...
        '_R' num2str(r) 'KM.'...
        fg_format];
    
    % Save figure  - output -
    folder_this = [folder_figs 'TIMESERIES/SSS/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name];
    
    % check fn existence
    fg_exist = exist(fg_name,'file');
    
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
end

%% [2.5.3] TIMESERIES: dSSS at each float location

% [dSSS = Argo - Satellite (NS and NM)]

% number of xtick
time_line = time_range(1):20:time_range(end);
time_line_diff = diff(time_line);

STATS_vars = {...
    ' mean dSSS_{NS}',' median dSSS_{NM}',...
    ' mean dSSS_{NM}',' median dSSS_{NM}',...
    ' mean dSSS_{NEMO}',' median dSSS_{NEMO}'};

NPLATS = unique(platform_irange);

for nn = 1: length(NPLATS)
    
    this_platform = NPLATS(nn);
    
    ind = find(platform_irange == this_platform);
    
    nprof = length(ind);

    X = time_range(ind);
    
    time_str1 = datestr(X(1),'yyyymmdd');
    time_str2 = datestr(X(end),'yyyymmdd');

    Y1_mean = dSSS_ns_irange_mean;
    Y2_mean = dSSS_nm_irange_mean;
    Y3_mean = dSSS_model_irange_mean;
    
    Y1_median = dSSS_ns_irange_median;
    Y2_median = dSSS_nm_irange_median;
    Y3_median = dSSS_model_irange_median;
    
    % =======================
    figure (25); clf; hold on
    
    h1 = plot(X,Y1_mean,'ks','MarkerFaceColor',color_lines(1,:),'MarkerSize',10); hold on
    h2 = plot(X,Y1_median,'-','color',color_lines(1,:),'LineWidth',5);
    h3 = plot(X,Y2_mean,'ko','MarkerFaceColor',color_lines(2,:),'MarkerSize',10); hold on
    h4 = plot(X,Y2_median,'-','color',color_lines(2,:),'LineWidth',5);
    h5 = plot(X,Y3_mean,'ko','MarkerFaceColor',color_lines(4,:),'MarkerSize',10); hold on
    h6 = plot(X,Y3_median,'-','color',color_lines(4,:),'LineWidth',5);
    
    
    plot(xlim, [0 0],'k--','linewidth',2)
    
    lg = legend([h1,h2,h3,h4,h5,h6],STATS_vars{:},...
        'fontsize',18,'location','SouthEast');
    
    ylim([-2 2])
    
    grid on
    grid minor
    box on
    
    set(gca,'FontSize',20)
    
    xticks(time_line)
    xticklabels(datestr(time_line,'mm'))
    % xtickangle(45)
    
    xlabel ('Time (months)','fontsize',20)
    ylabel ('dSSS [psu]','fontsize',20)
    
    title({[basin_str ' Salinity Argo (10 m) minus Satellite (dSSS) '];...
        ['ID: ' num2str(this_platform) ' (' num2str(nprof) ' profiles)'];
        [time_str1 '-' time_str2]})
    
    % Save figure  - output -
        fg_name = [...
            'ar_PLAT#' num2str(this_platform)...
            '_ARG2SAT_dSSS_TIMESERIES_'...
            time_str1 '_' time_str2...
            '_R' num2str(r) 'KM.'...
            fg_format];
    
    % Save figure  - output -
    folder_this = [folder_figs 'TIMESERIES/dSSS/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name];
    
    % check fn existence
    fg_exist = exist(fg_name,'file');
    
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
    
end


%% [2.5.4] TIMESERIES: ALL-dSSS at ALL-float location
% [dSSS = Argo - Satellite (NS and NM)]

% number of xtick
time_line = time_range(1):20:time_range(end);
time_line_diff = diff(time_line);

STATS_vars = {...
    ' mean dSSS_{NS}',' median dSSS_{NM}',...
    ' mean dSSS_{NM}',' median dSSS_{NM}',...
    ' mean dSSS_{NEMO}',' median dSSS_{NEMO}'};

X = time_range;

Y1_mean = dSSS_ns_irange_mean;
Y2_mean = dSSS_nm_irange_mean;
Y3_mean = dSSS_model_irange_mean;

Y1_median = dSSS_ns_irange_median;
Y2_median = dSSS_nm_irange_median;
Y3_median = dSSS_model_irange_median;


figure (25); clf; hold on

h1 = plot(X,Y1_mean,'ks','MarkerFaceColor',color_lines(1,:),'MarkerSize',10); hold on
h2 = plot(X,Y1_median,'-','color',color_lines(1,:),'LineWidth',5);
h3 = plot(X,Y2_mean,'ko','MarkerFaceColor',color_lines(2,:),'MarkerSize',10); hold on
h4 = plot(X,Y2_median,'-','color',color_lines(2,:),'LineWidth',5);
h5 = plot(X,Y3_mean,'ko','MarkerFaceColor',color_lines(4,:),'MarkerSize',10); hold on
h6 = plot(X,Y3_median,'-','color',color_lines(4,:),'LineWidth',5);


plot(xlim, [0 0],'k--','linewidth',2)

lg = legend([h1,h2,h3,h4,h5,h6],STATS_vars{:},...
    'fontsize',18,'location','SouthEast');

ylim([-2 2])

grid on
grid minor
box on

set(gca,'FontSize',20)

xticks(time_line)
xticklabels(datestr(time_line,'mm'))
% xtickangle(45)

xlabel ('Time (months)','fontsize',20)
ylabel ('dSSS [psu]','fontsize',20)

title({['Salinity Argo (10 m) minus Satellite '],...
    'Baltic ',...
    num2str(iyear)})

% Save figure  - output -
fg_name = ['ar_'...
    datestr(time_range(1),'yyyymmdd') '_'....
    datestr(time_range(end),'yyyymmdd') '_TIMESERIES_dSSS' ...
    '.'  fg_format];

% Save figure  - output -
folder_this = [folder_figs 'TIMESERIES/dSSS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];

% check fn existence
fg_exist = exist(fg_name,'file');

if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end


%% [2.6] TAYLOR DIAGRAMS

% Compute stats to include in Taylor diagram using Argo_10m 
% (averaged surface) as the reference measure 
% (use function: 'allstats')
% TAYLOR diagrram inputs:
%	[1] STDs: Standard deviations
%	[2] RMSs: Centered Root Mean Square Difference 
%	[3] CORs: Correlation

X_ref = SALT10_mean; % reference salinity

Y1 = sss_nm_irange_mn;
Y2 = sss_ns_irange_mn;
Y3 = sss_model_irange_mn;

% =========================
Y1_STATS = allstats(X_ref,Y1);
Y2_STATS = allstats(X_ref,Y2);
Y3_STATS = allstats(X_ref,Y3);

X_ref_STATS = Y1_STATS(:,1);

Y1_STATS = Y1_STATS(:,2);
Y2_STATS = Y2_STATS(:,2);
Y3_STATS = Y3_STATS(:,2);

% standard deviation
n = 2; 
STD_ALL = [...
    X_ref_STATS(n)...
    Y1_STATS(n)...
    Y2_STATS(n)...
    Y3_STATS(n)]; clear n

% Root Mean Square
n = 3; 
RMS_ALL = [...
    X_ref_STATS(n)...
    Y1_STATS(n)...
    Y2_STATS(n)...
    Y3_STATS(n)]; clear n

% Correlation
n = 4; 
CORR_ALL = [...
    X_ref_STATS(n)...
    Y1_STATS(n)...
    Y2_STATS(n)...
    Y3_STATS(n)]; clear n



% =================
figure
% output taylor
% hp: handle for each point
% ht: handle for each text
% axl: handle for axis
[hp ht axl] = taylordiag(STD_ALL,RMS_ALL,CORR_ALL);

marker_size = 20;

set(hp(1),'marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor','r','MarkerSize',marker_size)
set(hp(2),'marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color_lines(1,:),'MarkerSize',marker_size)
set(hp(3),'marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color_lines(2,:),'MarkerSize',marker_size)
set(hp(4),'marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color_lines(4,:),'MarkerSize',marker_size)


set(ht(1),'color','r','fontsize',18,'verticalalignment','bottom')
set(ht(2),'color','k','fontsize',18)
set(ht(3),'color','k','fontsize',18)
set(ht(4),'color','k','fontsize',18)

lg = legend(hp,{' [A] Argo_{10m} (reference)',' [B] NS',' [C] NM',' [D] NEMO'},...
    'Location','NorthEastOutside','fontsize',18);

% make legend square bigger 
% (https://uk.mathworks.com/matlabcentral/answers/420941-how-can-i-enlarge-the-height-of-legend-box)
lg.Position(1) = 0.80; % move legend left/right [0-1]
lg.Position(2) = 0.70; % move legend down/up [0-1] 
lg.Position(3) = 0.20; % move text in legend left/right [0-1] 
lg.Position(4) = 0.30; % increase/decrease legend box [0-1] 


% Save figure  - output -
fg_name = ['ar_'...
    datestr(time_range(1),'yyyymmdd') '_'....
    datestr(time_range(end),'yyyymmdd') '_TAYLOR_SSS' ...
    '.'  fg_format];

% Save figure  - output -
folder_this = [folder_figs 'TAYLOR/SSS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];

% check fn existence
fg_exist = exist(fg_name,'file');

if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end


% END-OF-SCRIPT





