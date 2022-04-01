%
%
% Syntax :
% >> as script: Baltic_SDN_BEC_binning_v1r0
%
% >> as function (not ready) [vars_out] = Baltic_SDN_BEC_binning(iyear,imonth,iregion,idata_type)
%
% Description
% Make bin (squares) of a given size around a reference (REF) data,
% to compute statistics (mean, mediam, std, number observations).
%
% Inputs
% iyear
% imonth
% flag_on (optional): if flag_on = 1 (default) it only keep data flagged
% as "good data" (flag = 49). Other wise (flag_on = 0), it does not apply
% flags.
%
% * for more informaiton about Baltic regions type "help Baltic_studyregion.m"
%
%
% Output
% vars_out is structure and Matlab (.mat) and NetCDF (.nc) files
% with binned data.
%
%
% current version: v1r0 (2020/07/15) -
%
% History
% v1r0 (2020/07/15) - original script to do colocations and stats.
%
% =========================================================================
% Author: rcatany
%
% =========================================================================

% clc;
close all

REF_str = 'SDN';
iyear = 2011:2013;
imonth = 9;
flag_on = 1;

% switch off warning (there is something about interpolation)
warning('off','all');

%% Flags script
run_test = 0;    % flag to run script (not to save figures);
save_output = 1; % flag to save output [1], or not [0];

bin_data = 1; % [bin_data = 1] bin REF, or not [bin_data = 0] - go to section 1.4

% bin spatial limits
ibin_sz = 25; % bind size (Km)
% convert bin size to degrees
bin_sz = ibin_sz/100; % *approx 1deg = 100 km, with high error at high lats

% bin temporal limmits
ndays = 9; % number of days contained in each BEC Argo file


if run_test == 1
    clc
    warning ('Baltic_seadatanet_BEC_colocation.m is in TEST MODE')
    
    disp(['RUNNING TEST MODE. Figures and files might NOT BE SAVED'])
    prompt_msm = 'Do you want to save script output? [1] yes, or [0] not ... ';
    save_output = input(prompt_msm);
end


%% flagging_data (glag = 49; SDN "good data"
% if isempty(varargin) || varargin{1} == 1
%     flag_on = 1;
%     flag_str = 'FLAG#49';
%     
% elseif varargin{1} == 0
%     flag_on = 0;
%     flag_str = 'FLAG#NO';
% else
%     flag_on = 1;
%     flag_str = 'FLAG#49';
%     
% end


%% set flag_str to notice filename
if flag_on == 1
    flag_str = 'FLAG#49';
elseif flag_on == 0
    flag_str = 'FLAG#00';
end


%% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance platform (in Km)
r_str = num2str(r);

% Areas at high-lats are about 25% bigger than at equator
lat_factor = 1.25;

% multiply lat_factor to get a more accurate sampling radius
r2 = r*lat_factor;

depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m

ibasin = 9; % Basin number: 9 (7 Arctic, 9 Baltic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/' REF_str '/']);

% ==================================================
%* Output variables to save
ref_vars_out = {...
    'SALT_raw','TEMP_raw','PRES_raw',...
    'lon_raw','lat_raw','time_raw',...
    'basin_str'};

ref_bin_out = {...
    'SALT_bin_mean','SALT_bin_median','SALT_bin_std',...
    'TEMP_bin_mean','TEMP_bin_median','TEMP_bin_std',...
    'PRES_bin',...
    'lon_bin','lat_bin','count_bin','time_bin','time_cnt','ibin_sz'};

save_vars_out = [ref_vars_out ref_bin_out];

% ==================================================
%% Loop each year Colocate ref-2-Satellite
% pre-locate vars

% pre-lcoate Stats dSSS = SATELLITE minus ARGO
nTOT    = 500;  % Total number of sat-to-argo colocations
ndepth  = 1;    % number of depth levels
grid_size = 1/4; % SMOS-grid size 1/4 (~0.25 km)
n = 9; % add extra elements to the matrix

% maximum number grid-points in irange + n extra elements
nele = (((r2/100)*4)/grid_size)+1+n;


%% Save output: binned product

% folder_out dependent of flag_on
% if flag_on == 0
%     folder_out = [path_root...
%         'SSS/' basin_str...
%         '/BEC/Validation/indata/' REF_str '/' REF_str '_mat' '/BINNED/monthly/'];
%     
% elseif flag_on == 1
%     folder_out = [path_root...
%         'SSS/' basin_str...
%         '/BEC/Validation/indata/' REF_str '/' REF_str '_mat' '/BINNED/' flag_str '/monthly/'];
% end

folder_out = [path_root...
    'SSS/' basin_str...
    '/BEC/Validation/indata/' REF_str '/' REF_str '_mat'...
    '/BINNED/' flag_str '/monthly/'];
foldercheck_raf(folder_out);

% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log REF_str '_MISSING_20200121.txt'];

folder_figs = '/Volumes/Rogue/scratch/Validation/';

folder_figs = [folder_figs basin_str '/' REF_str '/'];
if fg_save == 1
    foldercheck_raf(folder_figs); %! make folder_figs
end


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

SSS_QC = TT.SSS_QC;
SST_QC = TT.SST_QC;
PRES_QC = TT.depth_QC;
% =============================
%* Keep "SDN good data" only (flag = 49)
flag_good = 49;

if flag_on == 1
    ind = SSS_QC ~= flag_good & SST_QC ~= flag_good & PRES_QC ~= flag_good;
    SALT(ind) = NaN; TEMP(ind) = NaN;
end

platform_type = TT.metavar3; % [B]: Bottle; or [C]: CTD
platform_ID = TT.metavar5;
vars_measured = TT.metavar12;

z_grid = 0.5; % vertical interp at 0.5 m

% ===============================
%* function to make vertical intp (works with any prof., not only Argo)
[t_intp,s_intp,p_intp] = zinterp1_ARGO(TEMP,SALT,PRES,z_grid);

SALT = s_intp;
TEMP = t_intp;
PRES = p_intp;

ind = PRES > 100;

SALT(ind) = NaN;
TEMP(ind) = NaN;
PRES(ind) = NaN;

clear *intp

%% keep salinity values at depth_ref
depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m

nprof = length(lon);

SALT_ref = nan(1,nprof);
TEMP_ref = nan(1,nprof);

for ncount = 1:length(PRES(1,:))
    % Get REF salinity as value not deeper than depth_def (10 m). If
    % there are more than one measuremeent, take the median value
    SALT_ref_beta = SALT(find(PRES(:,ncount)<=depth_ref),ncount);
    TEMP_ref_beta = TEMP(find(PRES(:,ncount)<=depth_ref),ncount);
    
    ind = ~isnan(SALT_ref_beta);
    SALT_ref_beta = SALT_ref_beta(ind);
    TEMP_ref_beta = TEMP_ref_beta(ind);
    
    if length(SALT_ref_beta) > 1
        SALT_ref_beta = median(SALT_ref_beta);
        TEMP_ref_beta = median(TEMP_ref_beta);
        
    end
    
    if ~isnan(SALT_ref_beta)
        SALT_ref (ncount) = SALT_ref_beta;
        TEMP_ref (ncount) = TEMP_ref_beta;
    else
        S_beta = SALT(:,ncount);
        T_beta = TEMP(:,ncount);

        P_beta = PRES(:,ncount);
        
        ind = isnan(S_beta);
        S_beta(ind) = [];
        T_beta(ind) = [];
        P_beta(ind) = [];
        
        S_range = range(S_beta);
        T_range = range(T_beta);

        % keep profile: if [1] there are at least 2 salinity measures; [2]
        % within the first 30 m depth (we want to compare to surface); [3]
        % range bettween salintiies is lower than 0.1 (i.e. about satellite
        % accuracy)
        if length(S_beta) > 1 && S_range <= 0.1 && max(P_beta) < 30
            SALT_ref(ncount) = S_beta(end);
        else
            SALT_ref(ncount) = NaN;
        end
        
        % maxium temperature range (in the vertical) approximatelly 1
        % degree celcius
        if length(S_beta) > 1 && S_range <= 1 && max(P_beta) < 30
            TEMP_ref(ncount) = T_beta(end);
        else
            TEMP_ref(ncount) = NaN;
        end        
    end
    
end; clear *_beta*


%% [1.3] Grid SDN data to the smos-grid
%* ================================
%* ==== [NEW-CODE] : place in situ obs (REF) to the nearest satellite-node
%* store Baltic+ grid (lon/lat)
% folder_grid = ([path_root ...
%     'SSS/' basin_str '/BEC/' basin_str '_grid/']);
% 
% fn_grid = [folder_grid basin_str '_plusv1.0_grid.mat'];
% 
% %* load of a chosen satellite data (data_type) in a given date (itime_alfa)
% data_type = 5;
% itime_alfa = datenum(2012,09,01);
% 
% if exist(fn_grid,'file') ~= 2
%     [TT] = rd_smos_L4_BEC_v1r3(itime_alfa,ibasin,data_type);
%     
%     lon_sss = TT.lon;
%     lat_sss = TT.lat;
%     
%     % save Arctic grid
%     folder_grid = (['/Volumes/Rogue/Data/SSS/'...
%         basin_str '/BEC/' basin_str '_grid/']);
%     foldercheck_raf(folder_grid);
%     
%     save (fn_grid,'lon_sss','lat_sss');
% else
%     load (fn_grid)
% end


%* grid REF data to the satellite 'nearest' grid-node
% [lon_grid,lat_grid] = griddata_raf(lon,lat,lon_sss,lat_sss);

% lon = lon_grid;
% lat = lat_grid; clear *_grid

% =========
%* plot example of gridded and non-gridded data
plot_example = 0;

if plot_example == 1
    nn =1;
    run Baltic_seadatanet_BEC_colocations_PLOT.m
end

%* ===== [NEW] BIN REF data ============
%% [1.4] BIN reference data (space and time bins)

time_unique = unique(time_number);

for tt = 1:length(time_unique)
    
    % set binning temporal limits (ie. same as satellite Baltic+ data)
    time_cnt = time_unique(tt); % central time
    
    % define time limits in the bin (also filename)
    N = floor(ndays/2);
    time_start = time_cnt-N;
    time_end = time_cnt+N; clear N
    
    % time range in each binned file
    time_bin = time_start:time_end;
    
    time_start_str = datestr(time_start,'yyyymmdd');
    time_end_str = datestr(time_end,'yyyymmdd');
    
    % set file name out
    fn = [basin_str '_' REF_str ...
        '_BINNED_' num2str(ibin_sz) 'KM_'...
        time_start_str '_' time_end_str];
    
    fn_out = [folder_out fn '.mat'];
    
    % Do binning (if file does not exist)
    fn_exist = exist(fn_out,'file');
    
    if fn_exist ~= 2
        
        ind = find(time_number >= time_start & time_number<=time_end);
        
        lon_in = lon(ind);
        lat_in = lat(ind);
        
        SALT_in = SALT_ref(ind);
        TEMP_in = TEMP_ref(ind);
                        
        [TT1] = bindata_raf(ibasin,bin_sz,lon_in,lat_in,SALT_in);
        [TT2] = bindata_raf(ibasin,bin_sz,lon_in,lat_in,TEMP_in);
        
        lon_bin = TT1.lon_bin_center;
        lat_bin = TT1.lat_bin_center;
        
        count_bin = TT1.count_bin;
        
        SALT_bin_mean = TT1.mean_bin;
        SALT_bin_median = TT1.median_bin;
        SALT_bin_std = TT1.std_bin;
        
        TEMP_bin_mean = TT2.mean_bin;
        TEMP_bin_median = TT2.median_bin;
        TEMP_bin_std = TT2.std_bin;
        
        PRES_bin = ones(size(SALT_bin_mean)).*depth_ref;
        
        % ============
        % keep REF raw data
        SALT_raw = SALT_in;
        TEMP_raw = TEMP_in;
        PRES_raw = ones(size(SALT_in)).*depth_ref;
        lon_raw = lon_in;
        lat_raw = lat_in;
        time_raw = time_number(ind);
        
        time_cnt = ones(size(lon_raw)).*time_cnt;
        JULD = juliandate(datevec(time_cnt));
        
        % ============
        % save fn_out (binned ref files)
        save(fn_out,save_vars_out{:}); clear *_bin *_raw
        
    else
        
    end
end; clear tt

%% Save out in NetCDF format
% folder_cdf = [folder_out 'netCDF/'];
% foldercheck_raf(folder_out);
% 
% % Define filename
% fn_out_mat = fn_out; % matlab file
% fn_in = fn_out_mat;
% fn_exist = exist([fn_in '.mat'],'file');
% 
% if fn_exist == 2
%     %% Convert matlab file with binned data to netCDF    
%     fn_out = ([folder_out fn '.nc']);
%     
%     %% BUG writing NetCDF files
%     write_netcdf_SDN_Colocations_STATS(fn_in,fn_out)
%     
% end




%% funtion output

% TT.SALT_bin_mean = SALT_bin_mean;
% TT.SALT_bin_median = SALT_bin_median;
% TT.SALT_bin_std = SALT_bin_std;
% TT.TEMP_bin_mean = TEMP_bin_mean;
% TT.TEMP_bin_median = TEMP_bin_median;
% TT.TEMP_bin_std = TEMP_bin_std;
% 
% TT.coun_bin = count_bin;
% TT.lon_bin = lon_bin;
% TT.lat_bin = lat_bin;
% 
% vars = save_vars_out;
% N_vars = length(vars);
% 
% for x = 1:N_vars
% 
%     eval(['TT.' vars{x} ' = ' vars{x} ';'])
%     
% end
% 
% vars_out = TT;
