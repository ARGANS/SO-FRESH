
%
%
% Syntax :
% >> as script: Baltic_argo_BEC_binning_v1r0
%
% >> as function (not ready) [vars_out] = Baltic_argo_BEC_binning(iyear,imonth,iregion,idata_type)
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

clc;
close all;

REF_str = 'argo';
flag_on = 1; % bin data (flag = 1), or not (flag = 0);

% Enter time (year and months) to create binning data (ie. this can be
% vectors multiple years and months, or single years/months)
iYEARS = 2011:2013;
iMONTHS = 1:12;


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
    warning ('Baltic_argo_BEC_colocation.m is in TEST MODE')
    
    disp(['RUNNING TEST MODE. Figures and files might NOT BE SAVED'])
    prompt_msm = 'Do you want to save script output? [1] yes, or [0] not ... ';
    save_output = input(prompt_msm);
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
    'lon_raw','lat_raw','time_raw','platform',...
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
folder_out = [path_root...
    'SSS/' basin_str...
    '/BEC/Validation/indata/' REF_str '/argo_mat/BINNED/monthly/'];

foldercheck_raf(folder_out);

% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log REF_str '_MISSING_20200121.txt'];

folder_figs = '/Volumes/Rogue/scratch/Validation/';

folder_figs = [folder_figs basin_str '/' REF_str '/'];
if fg_save == 1
    foldercheck_raf(folder_figs); %! make folder_figs
end


%* plot example of gridded and non-gridded data
plot_example = 0;

if plot_example == 1
    nn =1;
    run Baltic_seadatanet_BEC_colocations_PLOT.m
end

%* ===== [NEW] BIN REF data ============
%% [1.4] BIN reference data (space and time bins)


% Define Year and Month
for yy = 1:length(iYEARS)
    iyear = iYEARS (yy);
    for mm = 1: length(iMONTHS)
        imonth = iMONTHS(mm);
        
        folder_in = ([folder_data 'argo_mat/' ...
            num2str(iyear) '/' sprintf('%02.0f',imonth) '/']);
        
        % ------------------------------------------------------------------------
        % DAYS IN A MONTH
        % In the next lines I am trying to make that script recognise the length of
        % the different folders. This is because each month has different number of
        % days. But I can not do it and I just do it manually because time is
        % running! (27/04/12)
        num_days = calendar(iyear,imonth)';
        num_days(num_days == 0) = [];
        num_days = length(num_days);
        
        % Days in a given month
        DAYS = 1:num_days;
        
        for dd = 1 : num_days
            iday = DAYS(dd);
            
            % set binning temporal limits (ie. same as satellite Baltic+ data)
            time_cnt = datenum(iyear,imonth,iday); % central time
            
            % =========================================
            % define time limits in the bin (also filename)
            N = floor(ndays/2);
            time_start = time_cnt-N;
            time_end = time_cnt+N; clear N
            % time range in each binned file
            time_bin = time_start:time_end;
            
            % Argo filename contains a time difference of 10-days instead
            % of 9-day difference (which is the time range set in satellite
            % data) -- BUG fix: add 1-day to filename IN (and inform BEC)
            time1_start_str = datestr(time_start,'yyyymmdd');
            time1_end_str = datestr(time_end+1,'yyyymmdd');
            
            time2_start_str = datestr(time_start,'yyyymmdd');
            time2_end_str = datestr(time_end,'yyyymmdd');
            
            
            
            %% Load Argo data (as provided by BEC (9-day files)
            % set filename IN
            fn1 = [REF_str '_' time1_start_str '_' time1_end_str];
            fn1_in = [folder_in fn1 '.mat'];
            
            % check fn_in existence
            fn1_exist = exist(fn1_in,'file');
            
            if fn1_exist == 2
                
                load(fn1_in)
                
                SALT = PSAL;
                PRES = ones(size(SALT)).*PRES;
                
                % number depth levels and number profiles in each file
                [ndepth,nprof] = size(SALT);
                
                time_number = datenum(YEAR,MONTH,DAY);
                % =========================================
                %% keep salinity values at depth_ref
                depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m
                for ncount = 1:nprof
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
                
                % =========================================
                % set file name OUT (file with binning data)
                fn2 = [basin_str '_' REF_str ...
                    '_BINNED_' num2str(ibin_sz) 'KM_'...
                    time2_start_str '_' time2_end_str];
                
                fn2_out = [folder_out fn2 '.mat'];
                
                % Do binning (if file does not exist)
                fn2_exist = exist(fn2_out,'file');
                
                
                
                if fn2_exist ~= 2
                    ind = find(time_number >= time_start & time_number<=time_end);
                    
                    if ~isempty(ind)
                        
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
                        JULD = juliandate(time_cnt);
                        
                        % ============
                        % save fn_out (binned ref files)
                        save(fn2_out,save_vars_out{:}); clear *_bin *_raw
                    end
                end
                
            end; clear TT
        end; clear dd % loop each day of the month
    end % for loop each MONTHS
end % for loop each YEARS

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
