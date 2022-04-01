% Syntax: Baltic_seadatanet_BEC_colocations.m (script)
%
% Description
% Make SeaDataNet-2-Satellite colocations whithin a searching
% distance and within
% a time window (i.e. ±7-days).
% [1] The output from this script is to feed into Baltic_argo_BEC_analyses.m
% [2] This script save colocations [Baltic_argo_BEC_colocations_VARS2SAVE.m]
%
% Input
% Seadatanet files prepared by BEC for validation SSS (Baltic+)
%
% Output
% Matlab files (.mat) with colocated profiles.
%
%
% current version: v1r0 (2020/02/21)
%
% =========================================================================
% Author: rcatany
%
% =========================================================================

clc;
close all

% switch off warning (there is something about interpolation)
warning('off','all');

%% Flags script
run_test = 0;    % flag to run script (not to save figures);
save_output = 1; % flag to save output [1], or not [0];

if run_test == 1
    clc
    warning ('Baltic_seadatanet_BEC_colocation.m is in TEST MODE')
    
    disp(['RUNNING TEST MODE. Figures and files might NOT BE SAVED'])
    prompt_msm = 'Do you want to save script output? [1] yes, or [0] not ... ';
    save_output = input(prompt_msm);
end

%% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance platform (in Km)

% Areas at high-lats are about 25% bigger than at equator
lat_factor = 1.25; 

% multiply lat_factor to get a more accurate sampling radius
r2 = r*lat_factor; 

depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m

fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/SeaDataNet/']);

%% Save Colocation output
folder_out = [path_root...
    'SSS/Baltic/BEC/Validation/indata/SeaDataNet/Colocations/'];
foldercheck_raf(folder_out);


% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'SDN_MISSING_20200121.txt'];

folder_figs = '/Volumes/Rogue/scratch/Validation/';

% ==========================
%
% iyear = 2011:2013;
% imonth = 1:12;


a = length(iyear);
b = length(imonth);

% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file

ibasin = 9; % Basin number: 7 (Arctic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

folder_figs = [folder_figs basin_str '/SeaDataNet/SDN2SAT/'];

if fg_save == 1
    foldercheck_raf(folder_figs); %! make folder_figs
end


% Output variables to save
SDN_vars_out = {...
    'lon_irange','lat_irange','time_irange',...
    'TEMP_intp','SALT_intp','PRES_intp',...
    'TEMP_irange','SALT_irange','PRES_irange',...
    'platform_type_irange','platform_ID_irange'};

sat_vars_out = {...
    'sss_nm_irange','sss_ns_irange','sss_gl_irange','sss_model_irange',...
    'sss_nm_error_irange','sss_ns_error_irange',...
    'lon_sss_irange','lat_sss_irange'};


colocated_vars_out = {...
    'sss_nm_irange_mn','sss_ns_irange_mn','sss_gl_irange_mn',...
    'sss_model_irange_mn',...
    'sss_nm_irange_md','sss_ns_irange_md','sss_gl_irange_md',...
    'sss_model_irange_md',...
    'sss_nm_error_irange_mn','sss_ns_error_irange_mn',...
    'sss_nm_error_irange_md','sss_ns_error_irange_md',...
    'sss_model_error_irange_mn','sss_model_error_irange_md'};

stats_X1_vars_out = {...
    'dSSS_nm_irange_nobs','dSSS_nm_irange_mean',...
    'dSSS_nm_irange_median','dSSS_nm_irange_std',...
    'dSSS_nm_irange_min','dSSS_nm_irange_max',...
    'dSSS_nm_irange_SEM','dSSS_nm_irange_Q'...
    'dSSS_nm_irange_Q2','dSSS_nm_irange_Q4',...
    'dSSS_nm_irange_IQR','dSSS_nm_irange_STATS_str'};

stats_X2_vars_out = {...
    'dSSS_ns_irange_nobs','dSSS_ns_irange_mean',...
    'dSSS_ns_irange_median','dSSS_ns_irange_std',...
    'dSSS_ns_irange_min','dSSS_ns_irange_max',...
    'dSSS_ns_irange_SEM','dSSS_ns_irange_Q'...
    'dSSS_ns_irange_Q2','dSSS_ns_irange_Q4',...
    'dSSS_ns_irange_IQR','dSSS_ns_irange_STATS_str'};

stats_X3_vars_out = {...
    'dSSS_gl_irange_nobs','dSSS_gl_irange_mean',...
    'dSSS_gl_irange_median','dSSS_gl_irange_std',...
    'dSSS_gl_irange_min','dSSS_gl_irange_max',...
    'dSSS_gl_irange_SEM','dSSS_gl_irange_Q'...
    'dSSS_gl_irange_Q2','dSSS_gl_irange_Q4',...
    'dSSS_gl_irange_IQR','dSSS_gl_irange_STATS_str'};


stats_X4_vars_out = {...
    'dSSS_model_irange_nobs','dSSS_model_irange_mean',...
    'dSSS_model_irange_median','dSSS_model_irange_std',...
    'dSSS_model_irange_min','dSSS_model_irange_max',...
    'dSSS_model_irange_SEM','dSSS_model_irange_Q'...
    'dSSS_model_irange_Q2','dSSS_model_irange_Q4',...
    'dSSS_model_irange_IQR','dSSS_model_irange_STATS_str'};

save_vars_out = [SDN_vars_out sat_vars_out colocated_vars_out...
    stats_X1_vars_out stats_X2_vars_out...
    stats_X3_vars_out stats_X4_vars_out];


%% Loop each year Colocate Argo-2-Satellite
% pre-locate vars

% pre-lcoate Stats dSSS = SATELLITE minus ARGO
nTOT    = 1000;  % Total number of sat-to-argo colocations
ndepth  = 200;   % number of depth levels

grid_size = 1/4; % SMOS-grid size 1/4 (~0.25 km)

n = 9; % add extra elements to the matrix

% maximum number grid-points in irange + n extra elements
nele = (((r2/100)*4)/grid_size)+1+n; 

ncount = 0; % set counter to zero


%% [1] Read SDN dataset

% Future releases might chnage filenanme of SDB dataset
fn_in = [folder_data...
    'data_from_SDN_2015-09_TS_BalticSea_QC_done_v2_filtered.nc'];


% function to read SDN dataset
[TT] = rd_SDN (fn_in);

lon = TT.lon;
lat = TT.lat;

time_number = TT.time_number;
time_str1 = datestr(time_number(1),'yyyymmdd');
time_str2 = datestr(time_number(end),'yyyymmdd');

SALT = TT.SSS_SDN;
TEMP = TT.SST_SDN;
PRES = TT.depth_SDN;

platform_type = TT.metavar3; % [B]: Bottle; or [C]: CTD
platform_ID = TT.metavar5;
vars_measured = TT.metavar12;

z_grid = 0.5; % vertical interp at 0.5 m

% function to make vertical intp (works with any prof., not only Argo)
[t_intp,s_intp,p_intp] = zinterp1_ARGO(TEMP,SALT,PRES,z_grid);

SALT = s_intp;
TEMP = t_intp;
PRES = p_intp;

ind = PRES > 100;

SALT(ind) = NaN;
TEMP(ind) = NaN;
PRES(ind) = NaN;

clear *intp


% ========================





%% Do colocations through time


for yy = 1:length(iyear)
    iYEAR = iyear(yy);
    
    fn_out = [folder_out 'SDN_'...
        sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
        sprintf('%02.0f',iYEAR)...
        '.mat'];
    
    % ==========================
    %% Pre-lcoate vars_out (ARGO-2-SAT colocation yearly files) -  snippet -
    run Baltic_SDN_BEC_colocations_VARS2SAVE.m % -  snippet -
    
    
    % =================================
    
    %% Colocate Argo-2-satellite (Loop through each month)
    fn_out_exist = exist(fn_out,'file');
    
    if fn_out_exist ~= 2 || run_test == 1
        
        disp(['Processing SDN-2-Satellite colocation '...
            basin_str ' ' num2str(iYEAR)]);
        

        
        
        for mm = 1:length(imonth)
            iMONTH = imonth(mm);
            
            disp(['year: ' num2str(iYEAR)...
                ' month: ' sprintf('%02.0f',iMONTH)])
            
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
                
                % >>>>>>>>>>>>>>>>>>>>
                % load SDN-matlab file [fn]
                %                 fn = (['argo_' ...
                %                     datestr(itime_start,'yyyymmdd') '_'...
                %                     datestr(itime_end,'yyyymmdd')]);
                %
                %                 fn_in = ([folder_in fn '.mat']);
                % >>>>>>>>>>>>>>>>>>>>
                
                % ===============
                %% [1] Load DATA
                %% [1.1] Load ARGO float
                if ismember(itime_start,time_number) == 1
                    
                    ind = time_number == itime_start;
                    
                    lon_alfa = lon(ind);
                    lat_alfa = lat(ind);
                    
                    platform_type_alfa = platform_type(ind);
                    platform_ID_alfa = platform_ID(ind);
                    
                    time_number_alfa = time_number(ind);
                    
                    SALT_alfa = SALT(:,ind);
                    TEMP_alfa = TEMP(:,ind);
                    PRES_alfa = PRES(:,ind);
                    
                    % station close to land are empty <??>
                    ind2 = any(SALT_alfa);
                    
                    lon_alfa(~ind2) = [];
                    lat_alfa(~ind2) = [];
                    time_number_alfa(~ind2) = [];
                    
                    SALT_alfa(:,~ind2) = [];
                    TEMP_alfa(:,~ind2) = [];
                    PRES_alfa(:,~ind2) = [];
                    
                    clear ind*
                    
                    if ~isempty(lon_alfa) && ~isempty(lat_alfa) && sum(any(SALT_alfa))~=0
                        % Cases with more than one float (profile) in one file
                        nprof1 = length(lon_alfa(:,1));
                        nprof2 = length(lat_alfa(:,1));
                        
                        % Argo lon-lat must be same size
                        if size(lon_alfa) == size(lat_alfa)
                            nprof = nprof1;
                        else
                            error(['Number of Argo longitudes and latitudes'...
                                ' must be equal'])
                        end
                        
                        % Make PRES to be same size as SALT (and TEMP)
                        [a1,b1] = size(SALT_alfa);
                        [a2,b2] = size(PRES_alfa);
                        
                        if b1 ~= b2
                            PRES_alfa = repmat(PRES_alfa,1,b1);
                        end
                        
                        %% [1.2] load Satellite SSS: Nodal Sampling (NS),
                        %% Nominal Sampling (NS) and Global product
                        
                        %% [1.2.1] Baltic-Nominal
                        data_type = 4;
                        [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                        
                        sss_nm = TT.sss;
                        sss_nm_error = TT.sss_error;
                        
                        %% [1.2.2] Baltic-Nodal Sampling
                        data_type = 5; % Baltic-Nodal Sampling
                        [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                        
                        sss_ns = TT.sss;
                        sss_ns_error = TT.sss_error;
                        
                        lon_sss = TT.lon;
                        lat_sss = TT.lat;
                        
                        %% [1.2.3] Baltic-Global product (v001)
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
                        
                        sss_gl = sss_gl2;
                        
                        clear sss_gl2 lon_gl lat_gl ind TT; % clear work space
                        
                        
                        %% [1.2.4] Load Model (NEMO)
                        
                        grid2baltic = 1; % grid model output to Baltic grid
                        
                        [YY,MM,DD] = datevec(itime_start);
                        
                        [TT] = nemo_rd_BEC_Baltic_FUN (YY,MM,DD,ibasin,grid2baltic);
                        
                        sss_model = TT.SSS_model;
                        lon_model = TT.lon;
                        lat_model = TT.lat;
                        
                        if lon_model ~= lon_sss
                            error('lon_model must be equal to lon_sss (satellite)')
                            
                        end
                        
                        if lat_model ~= lat_sss
                            error('lon_model must be equal to lon_sss (satellite)')
                       
                           
                        end                        
                        
                        % =================
                        
                        %% [1.3] Grid SDN data to the smos-grid
                        [lon_grid,lat_grid] = ...
                            griddata_raf(lon_alfa,lat_alfa,lon_sss,lat_sss);
                        
                        % plot example of gridded and non-gridded data
                        plot_example = 0;
                        
                        if plot_example == 1
                            nn =1;
                            run Baltic_seadatanet_BEC_colocations_PLOT.m
                        end
                        
                        
                        %% [2] Select Satellite-SSS data at Float location (within r distance in km)
                        for nn = 1: nprof
                            
                            % keep Argo lon/lat gridded to smos-grid
                            if nn == 1
                                lon_alfa = lon_grid;
                                lat_alfa = lat_grid;
                                clear lon_grid lat_grid lon_model lat_model
                            end
                            
                            DIST = Distance(lon_alfa(nn),lat_alfa(nn),lon_sss,lat_sss);
                            irange = find(DIST <= r2);
                            
                            Ln = length(irange);
                            
                            if Ln >= 1 
                                ncount = ncount+1;
                                
                                % to reset "repetition profile check"
                                rep_prof = -1;
                                
                                %% Satellite data at float location
                                sss_ns_irange (1:Ln,ncount) = sss_ns(irange);
                                sss_nm_irange (1:Ln,ncount) = sss_nm(irange);
                                sss_gl_irange (1:Ln,ncount) = sss_gl(irange);
                                
                                sss_ns_error_irange (1:Ln,ncount) = sss_ns_error(irange);
                                sss_nm_error_irange (1:Ln,ncount) = sss_nm_error(irange);
                                
                                sss_ns_irange_mn (ncount) = nanmean(sss_ns(irange));
                                sss_nm_irange_mn (ncount) = nanmean(sss_nm(irange));
                                sss_gl_irange_mn (ncount) = nanmean(sss_gl(irange));
                                
                                sss_ns_irange_md (ncount) = nanmedian(sss_ns(irange));
                                sss_nm_irange_md (ncount) = nanmedian(sss_nm(irange));
                                sss_gl_irange_md (ncount) = nanmedian(sss_gl(irange));
                                
                                sss_ns_error_irange_mn (ncount) = nanmean(sss_ns_error(irange));
                                sss_nm_error_irange_mn (ncount) = nanmean(sss_nm_error(irange));
                                
                                sss_ns_error_irange_md (ncount) = nanmedian(sss_ns_error(irange));
                                sss_nm_error_irange_md (ncount) = nanmedian(sss_nm_error(irange));
                                
                                lon_sss_irange (1:Ln,ncount) = lon_sss(irange);
                                lat_sss_irange (1:Ln,ncount) = lat_sss(irange);
                                
                                
                                %% Model (nemo) data at float location
                                sss_model_irange (1:Ln,ncount) = sss_model(irange);
                                sss_model_irange_mn (ncount) = nanmean(sss_model(irange));
                                sss_model_irange_md (ncount) = nanmedian(sss_model(irange));
                                
                                
                                %% [2.3] Argo Vertical profile interpolation
                                [a,b] = size(SALT_alfa);
                                
                                xROW = 1:a;
                                
                                % nprof_irange (1,ncount) = nprof;
                                
                                lon_irange (1,ncount) = lon_alfa(nn);
                                lat_irange (1,ncount) = lat_alfa(nn);
                                
                                time_irange (ncount) = time_number_alfa(nn);
                                
                                platform_type_irange(1,ncount) = char(platform_type_alfa(nn));
                                platform_ID_irange(1,ncount) = char(platform_ID_alfa(nn));
                                
                                PRES_irange (xROW,ncount) = PRES_alfa(:,nn);
                                SALT_irange (xROW,ncount) = SALT_alfa(:,nn);
                                TEMP_irange (xROW,ncount) = TEMP_alfa(:,nn);
                                
                                clear a b
                                
                                %% [2.3.1] Interpolate Z-direction (narrow depth levels)
                                % Interpolation specs: Grid dimensions
                                pres1 = 0;              % min pres bin
                                pres2 = max(PRES_alfa); % max depth (m)
                                z_grid = 0.5;           % interpolant distance (m) between depth levels
                                
                                TEMP1 = TEMP_alfa(:,nn);
                                SALT1 = SALT_alfa(:,nn) ;
                                PRES1 = PRES_alfa(:,nn) ;
                                
                                % interpolation function
                                [TEMP_beta_intp,SALT_beta_intp,PRES_beta_intp] =...
                                    zinterp1_ARGO(TEMP1,SALT1,PRES1,z_grid);
                                
                                [a,b] = size(TEMP_beta_intp);
                                
                                TEMP_intp(1:a,ncount) = TEMP_beta_intp;
                                SALT_intp(1:a,ncount) = SALT_beta_intp;
                                PRES_intp(1:a,ncount) = PRES_beta_intp;
                                
                                clear a b *_beta
                                
                                %% [3] dSSS = SDN - SATELLITE
                                % depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m
                                
                                SALT_10m = SALT_intp(find(PRES_intp(:,ncount)<=depth_ref),ncount);
                                
                                ind = ~isnan(SALT_10m);
                                SALT_10m = SALT_10m(ind);
                                
                                if length(SALT_10m) > 1
                                    SALT_10m = mean(SALT_10m);
                                end
                                
                                
                                if ~isempty (SALT_10m)
                                    sss_nm_irange_beta = sss_nm_irange(:,ncount);
                                    sss_ns_irange_beta = sss_ns_irange(:,ncount);
                                    sss_gl_irange_beta = sss_gl_irange(:,ncount);
                                    sss_model_irange_beta = sss_model_irange(:,ncount);
                                    
                                    nobs1 = sum(~isnan(sss_nm_irange_beta));
                                    nobs2 = sum(~isnan(sss_ns_irange_beta));
                                    nobs3 = sum(~isnan(sss_gl_irange_beta));
                                    nobs4 = sum(~isnan(sss_model_irange_beta));
                                    
                                    dSSS_nm_irange = sss_nm_irange_beta - SALT_10m;
                                    dSSS_ns_irange = sss_ns_irange_beta - SALT_10m;
                                    dSSS_gl_irange = sss_gl_irange_beta - SALT_10m;
                                    dSSS_model_irange = sss_model_irange_beta - SALT_10m;
                                    
                                    clear *_beta
                                    
                                end
                                
                                %% [3.1] Compute the Statistics of dSSS
                                
                                if ~isempty(SALT_10m)
                                    X1 = dSSS_nm_irange;
                                    X2 = dSSS_ns_irange;
                                    X3 = dSSS_gl_irange;
                                    X4 = dSSS_model_irange;
                                    
                                    % Number of Sat-2-Argo colocations
                                    dSSS_nm_irange_nobs(ncount) = nobs1;
                                    dSSS_ns_irange_nobs(ncount) = nobs2;
                                    dSSS_gl_irange_nobs(ncount) = nobs3;
                                    dSSS_model_irange_nobs(ncount) = nobs4;
                                    
                                    
                                    if any(~isnan(X1))
                                        [X1_mean,X1_median,...
                                            X1_std,X1_min,...
                                            X1_max,...
                                            X1_SEM,X1_Q,...
                                            X1_Q2,X1_Q4,...
                                            X1_IQR,X1_STATS_str] =...
                                            validation_stats(X1);
                                        
                                        % Keep STATS output
                                        % stats output 1- Nominal Sampling
                                        dSSS_nm_irange_mean(ncount) = X1_mean;
                                        dSSS_nm_irange_median(ncount) = X1_median;
                                        dSSS_nm_irange_std(ncount) = X1_std;
                                        dSSS_nm_irange_min(ncount) = X1_min;
                                        dSSS_nm_irange_max (ncount) = X1_max;
                                        dSSS_nm_irange_SEM (ncount) = X1_SEM;
                                        dSSS_nm_irange_Q(ncount,:) = X1_Q;
                                        dSSS_nm_irange_Q2(ncount) = X1_Q2;
                                        dSSS_nm_irange_Q4(ncount) = X1_Q4;
                                        dSSS_nm_irange_IQR(ncount) = X1_IQR;
                                        dSSS_nm_irange_STATS_str = X1_STATS_str;
                                        
                                    end
                                    
                                    
                                    if any(~isnan(X2))
                                        [X2_mean,X2_median,...
                                            X2_std,...
                                            X2_min,X2_max,...
                                            X2_SEM,X2_Q,...
                                            X2_Q2,X2_Q4,...
                                            X2_IQR,X2_STATS_str] =...
                                            validation_stats(X2);
                                        
                                        % stats output 2- Nodal Sampling
                                        dSSS_ns_irange_mean(ncount) = X2_mean;
                                        dSSS_ns_irange_median(ncount) = X2_median;
                                        dSSS_ns_irange_std(ncount) = X2_std;
                                        dSSS_ns_irange_min(ncount) = X2_min;
                                        dSSS_ns_irange_max (ncount) = X2_max;
                                        dSSS_ns_irange_SEM (ncount) = X2_SEM;
                                        dSSS_ns_irange_Q(ncount,:) = X2_Q;
                                        dSSS_ns_irange_Q2(ncount) = X2_Q2;
                                        dSSS_ns_irange_Q4(ncount) = X2_Q4;
                                        dSSS_ns_irange_IQR(ncount) = X2_IQR;
                                        dSSS_ns_irange_STATS_str = X2_STATS_str;
                                        
                                    end
                                    
                                    if any(~isnan(X3))
                                        [X3_mean,X3_median,...
                                            X3_std,...
                                            X3_min,X3_max,...
                                            X3_SEM,X3_Q,...
                                            X3_Q2,X3_Q4,...
                                            X3_IQR,X3_STATS_str] =...
                                            validation_stats(X3);
                                        % stats output 3- Global product
                                        dSSS_gl_irange_mean(ncount) = X3_mean;
                                        dSSS_gl_irange_median(ncount) = X3_median;
                                        dSSS_gl_irange_std(ncount) = X3_std;
                                        dSSS_gl_irange_min(ncount) = X3_min;
                                        dSSS_gl_irange_max (ncount) = X3_max;
                                        dSSS_gl_irange_SEM (ncount) = X3_SEM;
                                        dSSS_gl_irange_Q(ncount,:) = X3_Q;
                                        dSSS_gl_irange_Q2(ncount) = X3_Q2;
                                        dSSS_gl_irange_Q4(ncount) = X3_Q4;
                                        dSSS_gl_irange_IQR(ncount) = X3_IQR;
                                        dSSS_gl_irange_STATS_str = X3_STATS_str;
                                        
                                    end
                                    
                                    if any(~isnan(X4))
                                        [X4_mean,X4_median,...
                                            X4_std,...
                                            X4_min,X4_max,...
                                            X4_SEM,X4_Q,...
                                            X4_Q2,X4_Q4,...
                                            X4_IQR,X4_STATS_str] =...
                                            validation_stats(X4);
                                        % stats output 4- Model (NEMO) data
                                        dSSS_model_irange_mean(ncount) = X4_mean;
                                        dSSS_model_irange_median(ncount) = X4_median;
                                        dSSS_model_irange_std(ncount) = X4_std;
                                        dSSS_model_irange_min(ncount) = X4_min;
                                        dSSS_model_irange_max (ncount) = X4_max;
                                        dSSS_model_irange_SEM (ncount) = X4_SEM;
                                        dSSS_model_irange_Q(ncount,:) = X4_Q;
                                        dSSS_model_irange_Q2(ncount) = X4_Q2;
                                        dSSS_model_irange_Q4(ncount) = X4_Q4;
                                        dSSS_model_irange_IQR(ncount) = X4_IQR;
                                        dSSS_model_irange_STATS_str = X4_STATS_str;
                                        
                                    end
                                    
                                    clear X1_* X2_* X3_* SALT_10m
                                    
                                end
                            end
                        end; clear *_alfa *_beta*
                        
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
            
        end; clear mm dd TEMP SALT TEMP1 SALT1 lon lat
        
        
        %% Remove empty floats (*optimize space)
        
        % fixing bug: unwanted zeros in matrix
        ind = SALT_irange==0 & TEMP_irange==0; 
        
        SALT_irange(ind) = NaN;
        TEMP_irange(ind) = NaN; clear ind
        
        % REMEMBER save_vars_out:
        % save_vars_out = [SDN_vars_out colocated_vars_out...
        %                stats_X1_vars_out stats_X2_vars_out];
        
        [row1,col1] = size(SALT_intp);
        
%         if col1~= ncount
%             error('Number of columns of SALT_intp should be equal to ncounts')
%         else
%             nprof = ncount;
%         end
%         
        
        
        % index empty floats
        ind = find(all(isnan(SALT_irange),1) | isnan(time_irange)) ;
        
        if ~isempty (ind)
            % [BUG] Issue relating to the removing nan (columns/rows) --
            % fixed (badly) with a dirty for-loop and if-statements
            for nn = 1:length(save_vars_out)
                eval(['var_beta = ' save_vars_out{nn} ';']);
                
                % To remove platform profiles (NOT depth levels)
                [row2,col2] = size(var_beta);
                
                if col2 == col1 || col2 == nTOT
                    eval([save_vars_out{nn} '(:,ind) = [];']);
                    
                elseif row2 == col1 
                    eval([save_vars_out{nn} '(ind,:) = [];']);
                    
                elseif size(var_beta) == [nTOT,5]                    
                    eval([save_vars_out{nn} '(ind,:) = [];']);
                    
                end                
            end
        end 
        %clear row* col* ind
        
        %% save fn_out (with Argo-2-satellite colocations)
        disp({'New file with SDN-2-satellite colocations '; fn_out})
        save(fn_out,save_vars_out{:})
        
    elseif fn_out_exist == 2
        load(fn_out);
        
    else
        clc
        warning (['Missing Colocations file: ' fn_out])
        
        % write a log_file, recording the Argo (.mat) files
        msg_log = ([fn_in ' MISSING_FILE']);
        make_LOGFILE(fn_log,msg_log);
        
    end
end


% END-OF-SCRIPT





