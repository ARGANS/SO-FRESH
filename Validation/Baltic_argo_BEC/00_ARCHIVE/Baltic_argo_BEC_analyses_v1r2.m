% Syntax: Baltic_argo_BEC_analyses.m (script)
%
% Description
% Analyses of Argo floats as provided by BEC.
% Input: matlab files created from the original netCDF files (use
% ar_rd_BEC.m)
% 
% This script does the following:
% [1] Colocate Argo-2-satellite grid/time [Baltic_argo_colocations.m] 
% [2] Comparisons Argo Vs Satellite
% 
%
% Baltic Study [2011-2013]:
% *[1] BOXPLOTS SSS-Satellite Vs ARGO
% *[2] Argo Profiles (T and S) -- study column structure
% *[3]
%
% Data source:
% BEC products for validation
%
%
% current version: v1r2 (2020/01/20)
%
% History
% [3] Change order of script by:
%     --> 1/ Make daily plots Argo and SSS-satellite;
%     --> 2/ Save daily plots (to study TOP surface melting ice Vs. SSS)
%
% [2] Creation of snipet (functions) to make Argo vertical interpolations
% [1] creation of this script [20191114]
%
%
% =========================================================================

clc; clear;
close all

% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance platform (in Km)


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
iyear = 2013; % 2011:2012;
imonth = 7;   % 6:12;

a = length(iyear);
b = length(imonth);

% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file

ibasin = 9; % Basin number: 7 (Arctic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

folder_figs = [folder_figs basin_str '/Argo/Argo-BEC/'];

if fg_save == 1
    foldercheck_raf(folder_figs); %! make folder_figs
end

% pre-locate vars
Nobs_ALL = zeros(a,b);

Nobs_NA = zeros(a,b);
lon_NA = zeros(a,b);
lat_NA = zeros(a,b);
S_NA = zeros(a,b);
T_NA = zeros(a,b);

Nobs_PA = zeros(a,b);
lon_PA = zeros(a,b);
lat_PA = zeros(a,b);
S_PA = zeros(a,b);
T_PA = zeros(a,b);


% pre-lcoate Stats dSSS = SATELLITE minus ARGO
nTOT = 1000; % Total number of sat-to-argo colocations

TT = nan(nTOT,1);

% Nodal Sampling (X1)
X1_nobs = TT;
X1_mean = TT;
X1_median = TT;
X1_std = TT;
X1_min = TT;
X1_max = TT;
X1_SEM = TT;
X1_Q = nan(nTOT,5);
X1_Q2 = TT;
X1_Q4 = TT;
X1_IQR = TT;
X1_STATS_str = TT;

% Nominal Sampling (X2)
X2_nobs = TT;
X2_mean = TT;
X2_median = TT;
X2_std = TT;
X2_min = TT;
X2_max = TT;
X2_SEM = TT;
X2_Q = nan(nTOT,5);
X2_Q2 = TT;
X2_Q4 = TT;
X2_IQR = TT;
X2_STATS_str = TT;

ncount = 0; % set counter to zero

for yy = 1:length(iyear)
    for mm = 1:length(imonth)
        iYEAR = iyear(yy);
        iMONTH = imonth(mm);
        
        % folder name with BEC-ARGO data
        folder_in = ([folder_data...
            num2str(iYEAR) '/' sprintf('%02.0f',iMONTH) '/']);
        
        % number of days in month
        idays = calendar(iYEAR,iMONTH);
        idays (idays == 0) = [];
        idays = sort(idays);
        
        for dd = 1: length(idays)
            iDAY = idays(dd);
            itime_start = datenum(iYEAR,iMONTH,iDAY);
            itime_end = itime_start + ndays;
            
            iTIME = [itime_start itime_end]; % keep time in vector
            
            % load ARGO-matlab file [fn]
            fn = (['argo_' ...
                datestr(itime_start,'yyyymmdd') '_'...
                datestr(itime_end,'yyyymmdd')]);
            
            fn_in = ([folder_in fn '.mat']);
            
            % ===============
            %% [1] Load DATA
            %% [1.1] Load ARGO float
            if exist(fn_in,'file') == 2
                load (fn_in)
                
                if ~isempty(lon) && ~isempty(lat)
                    % Cases with more than one float (profile) in one file
                    nprof1 = length(lon(:,1));
                    nprof2 = length(lat(:,1));
                    
                    % Argo lon-lat must be same size
                    if size(lon) == size(lat)
                        nprof = nprof1;
                    else
                        error(['Number of Argo longitudes and latitudes'...
                            must be equal'])
                    end
                    
                    SALT = PSAL;
                    platform
                    
                    % Make PRES to be same size as SALT (and TEMP)
                    [a1,b1] = size(PSAL);
                    [a2,b2] = size(PRES);
                    
                    if b1 ~= b2
                        PRES = repmat(PRES,1,b1);
                    end
                    
                    %% [1.2] load Satellite SSS: Nodal Sampling (NS), Nominal Sampling (NS) and Global product 
                    
                    % 1/ Baltic-Nominal
                    data_type = 4;
                    [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                    
                    sss_nm = TT.sss;
                    sss_nm_error = TT.sss_error;
                    
                    % 2/ Baltic-Nodal Sampling
                    data_type = 5; % Baltic-Nodal Sampling
                    [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                    
                    sss_ns = TT.sss;
                    sss_ns_error = TT.sss_error;
                    
                    lon_sss = TT.lon;
                    lat_sss = TT.lat;
                    
                    % 3/ Baltic-Global product (v001)
                    data_type = 6;
                    [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                    
                    sss_gl = TT.sss;
                    
                    lon_gl = TT.lon;
                    lat_gl = TT.lat;
                    
                    % homogenize grid (global product to baltic+)
                    [sss_gl2] = griddata(lon_gl,lat_gl,sss_gl,lon_sss,lat_sss);
                    
                    ind = lon_sss >= xmin & lon_sss <= xmax |...
                        lat_sss >= ymin & lat_sss <= ymax ;
                    
                    sss_gl2 (ind == 0) = NaN;
                    
                    sss_gl = sss_gl2; clear sss_gl2 lon_gl lat_gl ind; % clear work space
                    
                    % =================
                    
                    for nn = 1: nprof
                        ncount = ncount+nn;
                        %% [2] Select Satellite-SSS data at Float location (within r distance in km)
                        % r = 25;
                        
                        DIST = Distance(lon(nn),lat(nn),lon_sss,lat_sss);
                        irange = find(DIST <= r);
                        
                        Ln = length(irange);
                        
                        if Ln >= 1
                            sss_ns_irange = sss_ns(irange);
                            sss_nm_irange = sss_nm(irange);
                            sss_gl_irange = sss_gl(irange);
                            
                            sss_ns_error_irange = sss_ns_error(irange);
                            sss_nm_error_irange = sss_nm_error(irange);
                            
                            sss_ns_irange_mn = nanmean(sss_ns(irange));
                            sss_nm_irange_mn = nanmean(sss_nm(irange));
                            sss_gl_irange_mn = nanmean(sss_gl(irange));
                            
                            sss_ns_irange_md = nanmedian(sss_ns(irange));
                            sss_nm_irange_md = nanmedian(sss_nm(irange));
                            sss_gl_irange_md = nanmedian(sss_gl(irange));
                            
                            sss_ns_error_irange_mn = nanmean(sss_ns_error(irange));
                            sss_nm_error_irange_mn = nanmean(sss_nm_error(irange));
                            
                            sss_ns_error_irange_md = nanmedian(sss_ns_error(irange));
                            sss_nm_error_irange_md = nanmedian(sss_nm_error(irange));
                                         
                            lon_sss_irange = lon_sss(irange);
                            lat_sss_irange = lat_sss(irange);
                            
                            % ===============
                            %% [2.1] boxplots: Satellite data at float location
                            % shall we use the median of SSS to be compared againts
                            
                            % string product name
                            PRODS = {'SSS_{NS}','SSS_{NM}','SSS_{gl}'};
                            
                            
                            make_plot = 1; % flag to make plot [1], or not [0]
                            if make_plot == 1
                            
                            figure(1); clf; hold on;
                            set(gcf,'DefaultAxesFontSize',24);
                            set(gca,'TickLabelInterpreter','tex');
                            
                            bplot = boxplot...
                                ([sss_ns_irange,sss_nm_irange,sss_gl_irange],...
                                'labels',PRODS);
                            
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
                            
                            end
                            
                            %% *[2.2] Plot: Argo geolocation
                            % Plots in one map locations of same platform
                            
                            make_plot = 0; % flag to make plot [1], or not [0]
                            if make_plot == 1

                            % snippet to plot map and zoom-in to Argo location
                            plot_map_argo_BEC
                            
                            end
                            
                            %% *[2.3] Plot Argo profiles
                            
                            % Interpolate Z-direction (narrow depth levels)
                            % Interpolation specs: Grid dimensions
                            pres1 = 0;              % min pres bin
                            pres2 = max(PRES);      % max depth (m)
                            z_grid = 0.5;           % interpolant distance (m) between depth levels
                            
                            TEMP1 = TEMP (:,nn);
                            SALT1 = PSAL (:,nn);
                            PRES1 = PRES (:,nn);
                            
                            [t_intp,s_intp,p_intp] =...
                                zinterp1_ARGO(TEMP1,SALT1,PRES1,z_grid);
                            
                            %% Example
                            %% [2.3.1] plot raw Vs interp T and S of a given Profile (nprof)
                            TEMP2 = t_intp(:,1);
                            SALT2 = s_intp(:,1);
                            PRES2 = p_intp(:,1);
                            
                            str1 = 'RAW';
                            str2 = 'intp';
                            
                            
                            make_plot = 0; % flag to make plot [1], or not [0]
                            if make_plot == 1
                                
                            figure(4); clf; hold on;
                            set(gcf,'DefaultAxesFontSize',18);
                            
                            plot_profile_argo_BEC(TEMP1,SALT1,PRES1,...
                                TEMP2,SALT2,PRES2,platform(nn),...
                                iTIME,str1,str2)
                            
                            % Save figure  - output -
                            fg_name = [fn '_RAWvsINTP_PLAT'...
                                num2str(platform(nn))...
                                '_PROF#' num2str(ID(nn)) ...
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
                            
                            %% [2.3.2] More than one profile-platform, make mean and compare all profiles
                            if nprof > 1
                                [t_intp,s_intp,p_intp] =...
                                    zinterp1_ARGO(TEMP,SALT,PRES,z_grid);
                                
                                TEMP1 = t_intp;
                                SALT1 = s_intp;
                                PRES1 = p_intp;
                                
                                TEMP2 = nanmean(t_intp,2);
                                SALT2 = nanmean(s_intp,2);
                                PRES2 = p_intp(:,1);
                                
                                str1 = 'all-profiles';
                                str2 = 'mean';
                                
                                % ===========
                                
                                make_plot = 0;
                                
                                if make_plot == 1
                                figure(4); clf; hold on
                                set(gcf,'DefaultAxesFontSize',18);
                                
                                plot_profile_argo_BEC(TEMP1,SALT1,PRES1,...
                                    TEMP2,SALT2,PRES2,platform(nn),...
                                    iTIME,str1,str2)
                                
                                % Save figure  - output -
                                fg_name = [fn '_BOUNDEDPLOT_PLAT'...
                                    num2str(platform(nn)) ...
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
                            
                            end
                            
                            
                            %% [3] dSSS = ARGO - SATELLITE
                            ind = find (PRES2>= 8 & PRES2 <=12);
                            
                            SALT_10m = SALT2(find(PRES2>=10,1,'first'));
                            
                            nobs1 = sum(~isnan(sss_ns_irange));
                            nobs2 = sum(~isnan(sss_nm_irange));
                            
                            dSSS_ns_irange = sss_ns_irange - SALT_10m;
                            dSSS_nm_irange = sss_nm_irange - SALT_10m;
                            
                            %% [3.1] Compute the Statistics of dSSS

                            X1 = dSSS_ns_irange;
                            X2 = dSSS_nm_irange;
                            
                            % Number of Sat-2-Argo colocations 
                            X1_nobs(ncount) = nobs1;
                            X2_nobs(ncount) = nobs2;
                                                        
                            [X1_mean(ncount,1),X1_median(ncount,1),...
                                X1_std(ncount,1),X1_min(ncount,1),X1_max(ncount,1),...
                                X1_SEM(ncount,1),X1_Q(ncount,:),X1_Q2(ncount),X1_Q4(ncount,1),...
                                X1_IQR(ncount,1),X1_STATS_str] =...
                                validation_stats(X1);
                            
                            [X2_mean(ncount),X2_median(ncount),X2_std(ncount),...
                                X2_min(ncount),X2_max(ncount),...
                                X2_SEM(ncount),X2_Q(ncount,:),X2_Q2(ncount),X2_Q4(ncount),...
                                X2_IQR(ncount),X2_STATS_str] =...
                                validation_stats(X2);
%                            
                        end
                    end
                    
                elseif isempty(lon) || isempty(lat)
                    msg_log = ([fn_in ' EMPTY_FILE']);
                    make_LOGFILE(fn_log,msg_log);
                end
                
            elseif exist(fn_in,'file') == 0
                % write a log_file, recording the Argo (.mat) files
                msg_log = ([fn_in ' MISSING_FILE']);
                make_LOGFILE(fn_log,msg_log);
            end
        end
    end
end


%%
ind1 = isnan(X1_nobs);
ind2 = isnan(X2_nobs);

% Nodal Sampling (X1)
X1_nobs (ind1) = [];
X1_mean (ind1) = [];
X1_median (ind1) = [];
X1_std (ind1) = [];
X1_min (ind1) = [];
X1_max (ind1) = [];
X1_SEM (ind1) = [];
X1_Q (ind1,:) = [];
X1_Q2 (ind1) = [];
X1_Q4 (ind1) = [];
X1_IQR (ind1) = [];
X1_STATS_str;

% Nominal Sampling (X2)
X2_nobs (ind2) = [];
X2_mean (ind2) = [];
X2_median (ind2) = [];
X2_std (ind2) = [];
X2_min (ind2) = [];
X2_max (ind2) = [];
X2_SEM (ind2) = [];
X2_Q (ind2,:) = [];
X2_Q2 (ind2) = [];
X2_Q4 (ind2) = [];
X2_IQR (ind2) = [];
X2_STATS_str;

%%
STATS_vars = {'mean dSSS_{NS}','median dSSS_{NM}',...
    'mean dSSS_{NM}','median dSSS_{NM}'};

figure

h1 = plot(X1_mean,'ks','MarkerFaceColor','b','MarkerSize',10); hold on
h2 = plot(X1_median,'b-','LineWidth',5);
h3 = plot(X2_mean,'ko','MarkerFaceColor','r','MarkerSize',10); hold on
h4 = plot(X2_median,'r-','LineWidth',5);

plot(xlim, [0 0],'k--','linewidth',2)

lg = legend([h1,h2,h3,h4],STATS_vars{:},'fontsize',18,'location','SouthEast');

ylim([-2 0.5])

grid on
grid minor
box on

set(gca,'FontSize',20)
xlabel ('Colocation numbers','fontsize',20)
ylabel ('dSSS [psu]','fontsize',20)

title({['Salinity Argo (10 m) minus Satellite '],...
    'Baltic ',... 
    num2str(iyear)})

% Save figure  - output -
fg_name = [fn '_BOUNDEDPLOT_dSSS' ...
    '.'  fg_format];

% Save figure  - output -
folder_this = [folder_figs 'dSSS/'];

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





