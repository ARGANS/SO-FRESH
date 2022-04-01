function [vars_out] = Baltic_argo_BEC_colocations_v2r0(iyear,imonth,iregion,idata_type)
%
% Syntax (function):
% [vars_out] = Baltic_argo_BEC_colocations(iyear,imonth,iregion,idata_type)
%
% Description
% Make Argo-2-Satellite colocations whithin a searching distance and within
% a time window (i.e. ±7-days).
% [1] The output from this script is to feed into Arctic_argo_BEC_analyses.m
% [2] This script save colocations [Arctic_argo_BEC_colocations_VARS2SAVE.m]
%
% Argo files were prepared by BEC for validation SSS (Baltic+ and Arctic+)
%
%
% Input
% iyear, imonth
% iregion: There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
% [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
% [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
% [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
% [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)
% [5] ALL-Baltic             [ALL] 
% 
% * for more informaiton about Baltic regions type "help Baltic_studyregion.m"
% 
%
%
% data_type:   [1] Nominal, [3] Nodal Sampling,
%              [6] Global BEC (v001);
%              [7] CCI+SSS (v01.07);
%              [8] Model (NEMO) data in Baltic
%
%
% Output
% vars_out is structure and Matlab files (.mat) with colocated profiles.
%
%
% current version: v2r0 (2020/03/18) -
%         [1] converted script into function;
%         [2] split up colocations and stats in
%             two different scripts. This script does colocation (only)
%             to input in stats.
%
% History
% v1r2 (2020/01/20) - original script to do colocations and stats
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
    warning ([basin_str '_argo_BEC_colocation.m is in TEST MODE'])
    
    disp(['RUNNING TEST MODE. Figures and files might NOT BE SAVED'])
    prompt_msm = 'Do you want to save script output? [1] yes, or [0] not ... ';
    save_output = input(prompt_msm);
end


% set the data_type
if idata_type == 1
    version_str = 'plusv3r0';
elseif idata_type == 3
    version_str = 'v2r0';
elseif idata_type == 4
    version_str = 'v1r0_NOM';
elseif idata_type == 5
    version_str = 'v1r0_NS';
elseif idata_type == 6
    version_str = 'v001';
elseif idata_type == 7
    version_str = 'CCI+SSSv1.7';
elseif idata_type == 8
    version_str = 'REANALYSIS_PHY_003_011';
end

%% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance platform (in Km)
lat_factor = 1.25; % Areas at high-lats are about 25% bigger than at equator

r_str = num2str(r);

r2 = r*lat_factor; % multiply lat_factor to get a more accurate sampling radius

depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m

ibasin = 9; % Basin number: 9 (7 Arctic, 9 Baltic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);


path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/' basin_str '/BEC/Validation/indata/argo/argo_mat/']);

%% Save Colocation output
folder_out = [path_root...
    'SSS/' basin_str...
    '/BEC/Validation/indata/argo/Colocations/monthly/' version_str '/'];
foldercheck_raf(folder_out);

%% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = ['/Volumes/Rogue/Data/SSS/' basin_str '/BEC/Validation/indata/'];
fn_log = [folder_log 'ARGO_MISSING_20200121.txt'];

% ==========================

% get dimension of var_out (a,b)
a = length(iyear);
b = length(imonth);

% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file


% Output variables to save
argo_vars_out = {...
    'platform_irange','ID_irange',...
    'PRES_irange','SALT_irange','TEMP_irange',...
    'lon_irange','lat_irange','JULD_irange','time_irange'...
    'TEMP_intp','SALT_intp','PRES_intp','region','basin_str'};

sat_vars_out = {...
    'sss_irange',...
    'sss_error_irange',...
    'lon_sss_irange','lat_sss_irange','version_str'};

colocated_vars_out = {...
    'sss_irange_mn',...
    'sss_irange_md',...
    'sss_error_irange_mn'};

save_vars_out = [argo_vars_out sat_vars_out colocated_vars_out];

%% Loop each year Colocate Argo-2-Satellite
% pre-locate vars

% pre-lcoate Stats dSSS = SATELLITE minus ARGO
nTOT    = 500;   % Total number of sat-to-argo colocations
ndepth  = 200;   % number of depth levels

grid_size = 1/4; % SMOS-grid size 1/4 (~0.25 km)

n = 9; % add extra elements to the matrix
nele = (((r2/100)*4)/grid_size)+1+n; % maximum number grid-points in irange + n extra elements

ncount = 0; % set counter to zer

% Colocatiion fn_out for each Regional study in the Baltic
if strcmpi(iregion,'ARB') || iregion == 1% Arkona Basin
    region = 'ARB';
    
elseif strcmpi(iregion,'BOS')|| iregion == 2 % Bothian Sea
    region = 'BOS';
    
elseif strcmpi(iregion,'GOF')|| iregion == 3 % Gulf of Finland
    region = 'GOF';
    
elseif strcmpi(iregion,'NBP')|| iregion == 4 % Gulf of Finland
    region = 'NBP';
    
elseif strcmpi(iregion,'ALL') || iregion == 5 % all-Baltic region stats
    region = 'ALL';
    
end


% loop through time (year, month and days) to load and colocate argo

for yy = 1:length(iyear)
    iYEAR = iyear(yy);
    
    disp(['Processing Argo-2-Satellite colocation '...
        basin_str ' ' num2str(iYEAR)]);
    
    for mm = 1:length(imonth)
        iMONTH = imonth(mm);
        
        % =============================
        % fn_out (.mat)
        fn_out = [folder_out basin_str '_' region '_argo_R' r_str '_'...
            sprintf('%02.0f',depth_ref) 'm_COLOCATIONS_'...
            sprintf('%02.0f',iYEAR) sprintf('%02.0f',iMONTH) '_'...
            version_str '.mat'];
        
        %% Colocate Argo-2-satellite (Loop through each month)
        fn_out_exist = exist(fn_out,'file');
        
        if fn_out_exist ~= 2 || run_test == 1
            
            disp(['year: ' num2str(iYEAR)...
                ' month: ' sprintf('%02.0f',iMONTH)])
            
            % ==========================
            %% Pre-lcoate vars_out (ARGO-2-SAT colocation yearly files) -  snippet -
            % run Baltic_argo_BEC_colocations_VARS2SAVE.m % -  snippet -
            run ([basin_str '_argo_BEC_colocations_VARS2SAVE.m']) % -  snippet -
            
            % =================================
            % folder name with BEC-ARGO data
            folder_in = ([folder_data...
                num2str(iYEAR) '/' sprintf('%02.0f',iMONTH) '/']);
            
            % number of days in month
            idays = calendar(iYEAR,iMONTH);
            idays (idays == 0) = [];
            idays = sort(idays);
            
            for dd = 1:length(idays)
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
                    
                    if exist('lon','var') && exist('lat','var')
                        
                        % Regional study in the Baltic
                        if strcmpi(region,'ARB') % North Atlantic
                            
                            ind_reg = lon >=15 & lon <= 18;
                            
                        elseif strcmpi(region,'BOS')
                            
                            ind_reg = lon>=18 & lon <= 23 & lat > 60 & lat <= 64;
                            
                        elseif strcmpi(region,'GOF')
                            
                            ind_reg = lon>=23 & lon <= 26 & lat > 58 & lat <= 62;
                            
                        elseif strcmpi(region,'NBP')
                            
                            ind_reg = lon>=18 & lon <= 23 & lat > 56 & lat <= 59;
                            
                        elseif strcmpi(region,'ALL') % all-Arctic region stats
                            ind_reg = ones(size(lon));
                            
                        end
                        
                        lon(~ind_reg) = [];
                        lat(~ind_reg) = []; clear ind_reg
                        
                        % ===================================
                        if ~isempty(lon) && ~isempty(lat)
                            % Cases with more than one float (profile) in one file
                            nprof1 = length(lon(:,1));
                            nprof2 = length(lat(:,1));
                            
                            % Argo lon-lat must be same size
                            if size(lon) == size(lat)
                                nprof = nprof1;
                            else
                                error(['Number of Argo longitudes and latitudes'...
                                    ' must be equal'])
                            end
                            
                            SALT = PSAL;
                            platform;
                            
                            % Make PRES to be same size as SALT (and TEMP)
                            [a1,b1] = size(PSAL);
                            [a2,b2] = size(PRES);
                            
                            if b1 ~= b2
                                PRES = repmat(PRES,1,b1);
                            end
                            
                            %% [1] load Satellite SSS:
                            %% [1] Baltic+ Nominal, [2] Baltic Nodal Sampling, [3] Global BEC, [4] CCI+SSS
                            %% [1.2.1] Baltic+ [BEC-SSS] Nominal v1.0 20110101 - 20131227
                            
                            if idata_type == 4 && itime_start <= datenum(2013,12,27)
                                data_type = idata_type; % Baltic+ (SSS-BEC NM v1.0)
                                
                                [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                                
                                lon_sss1 = TT.lon;
                                lat_sss1 = TT.lat;
                                
                                sss = TT.sss;
                                sss_error = TT.sss_error;
                                
                                lon_sss = lon_sss1;
                                lat_sss = lat_sss1;
                                
                            elseif idata_type == 1 && itime_start > datenum(2019,08,31)
                                sss = nan;
                                lon_sss = NaN;
                                lat_sss = NaN;
                            end
                            
                            
                            
                            
                            %% [1.2.2] Baltic+ Nodal Sampling [version v1.0: 20110101 - 20180102]
                            
                            if idata_type == 5 && itime_start <= datenum(2013,12,27)
                                data_type = idata_type; % Baltic+ (SSS-BEC NS v1.0)
                                
                                [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                                
                                lon_sss1 = TT.lon;
                                lat_sss1 = TT.lat;
                                
                                sss = TT.sss;
                                sss_error = TT.sss_error;
                                
                                lon_sss = lon_sss1;
                                lat_sss = lat_sss1;
                                
                            elseif idata_type == 1 && itime_start > datenum(2019,08,31)
                                sss = nan;
                                lon_sss = NaN;
                                lat_sss = NaN;
                            end
                            
                            
                            
                            %% load Baltic+ grid -to convert all products to the same grid
                            if idata_type == 6 || idata_type == 7 || idata_type == 8
                                
                                data_type = 5;
                                % store Baltic+ grid (lon/lat)
                                folder_grid = ([path_root ...
                                    'SSS/' basin_str '/BEC/Baltic_grid/']);
                                
                                fn_grid = [folder_grid 'Baltic_plusv1.0_grid.mat'];
                                
                                if exist(fn_grid,'file') ~= 2
                                    [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                                    
                                    lon_sss = TT.lon;
                                    lat_sss = TT.lat;
                                    
                                    % save Arctic grid
                                    folder_grid = (['/Volumes/Rogue/Data/SSS/Baltic/BEC/' basin_str '_grid/']);
                                    foldercheck_raf(folder_grid);
                                    
                                    save (fn_grid,'lon_sss','lat_sss');
                                    
                                else
                                    load (fn_grid)
                                end
                            end
                            
                            %% [1.2.3] Baltic-Global product (v001)
                            
                            if idata_type == 6 && iYEAR < 2010
                                data_type = 6;
                                [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                                
                                sss_beta = TT.sss;
                                
                                lon_beta = TT.lon;
                                lat_beta = TT.lat;
                                
                                % homogenize grid (global product to Baltic+)
                                [sss_beta2] = griddata(lon_beta,lat_beta,sss_beta,lon_sss,lat_sss);
                                
                                ind = lon_sss >= xmin & lon_sss <= xmax |...
                                    lat_sss >= ymin & lat_sss <= ymax ;
                                
                                sss_beta2 (ind == 0) = NaN;
                                
                                sss = sss_beta2;
                                sss_error = nan(size(sss));
                                
                                clear *_beta* ind TT; % clear work space
                                
                            elseif idata_type == 6 && iYEAR >=2010
                                sss = nan;
                                lon_sss = NaN;
                                lat_sss = NaN;
                                
                                
                            end
                            
                            
                            %% [1.2.4] CCI+SSS data (v01.07)
                            % CCI+SSS data: 20100106 - 20181101
                            if idata_type == 7 && itime_start <= datenum(2018,11,01)
                                
                                plot_cci_ex = 0;
                                cci_product = 1; % [1] '7-days', or [2] '30-days' product
                                [TT] = rd_sss_cci(itime_start,ibasin,cci_product,plot_cci_ex);
                                
                                sss_beta = TT.sss;
                                
                                sss_error_beta = TT.sss_bias_std;
                                
                                % get lat/lon from Baltic+
                                lon_sss1 = lon_sss;
                                lat_sss1 = lat_sss;
                                
                                lon_beta = TT.lon;
                                lat_beta = TT.lat;
                                
                                % Grid Baltic product to the same grid (all in
                                % Baltic+ grid)
                                sss = griddata(lon_beta,lat_beta,sss_beta,lon_sss1,lat_sss1);
                                sss_error = griddata(lon_beta,lat_beta,sss_error_beta,lon_sss1,lat_sss1);
                                
                                lon_sss = lon_sss1;
                                lat_sss = lat_sss1;
                                
                                clear lon_sss1 lat_sss1 *_beta*
                                
                            elseif idata_type == 7 && itime_start > datenum(2018,11,01)
                                sss = nan;
                                lon_sss = NaN;
                                lat_sss = NaN;
                                
                            end
                            
                            
                            %% [1.2.5] Load Model in Baltic
                            
                            if idata_type == 8 && itime_start <= datenum(2013,12,27)
                                
                                grid2baltic = 1; % grid model output to Baltic grid
                                [YY,MM,DD] = datevec(itime_start);
                                
                                
                                data_type = idata_type; % Baltic+ MODEL
                                
                                [TT] = nemo_rd_BEC_Baltic_FUN (YY,MM,DD,ibasin,grid2baltic);
                                
                                sss_beta = TT.SSS_model;
                                
                                % get lat/lon from Baltic+
                                lon_sss1 = lon_sss;
                                lat_sss1 = lat_sss;
                                
                                lon_beta = TT.lon;
                                lat_beta = TT.lat;
                                
                                % Grid Baltic product to the same grid (all in
                                % Baltic+ grid)
                                sss = griddata(lon_beta,lat_beta,sss_beta,lon_sss1,lat_sss1);
                                
                                sss_error = nan(size(sss));
                                
                                lon_sss = lon_sss1;
                                lat_sss = lat_sss1;
                                
                                clear lon_sss1 lat_sss1 *_beta*
                                
                            elseif idata_type == 8 && itime_start > datenum(2013,12,27)
                                sss = nan;
                                lon_sss = NaN;
                                lat_sss = NaN;
                            end
                            
                            
                            % =================
                            
                            %% [1.3] Grid Argo data to the smos-grid
                            [lon_grid,lat_grid] = ...
                                griddata_raf(lon,lat,lon_sss,lat_sss);
                            
                            % plot example of gridded and non-gridded data
                            plot_example = 0;
                            
                            if plot_example == 1
                                nn =1;
                                run Baltic_argo_BEC_colocations_PLOT_griddedARGO.m
                            end
                            
                            
                            %% [2] Select Satellite-SSS data at Float location (within r distance in km)
                            for nn = 1: nprof
                                
                                % keep Argo lon/lat gridded to smos-grid
                                if nn == 1
                                    lon = lon_grid;
                                    lat = lat_grid;
                                    clear lon_grid lat_grid
                                end
                                
                                % Repetion profile - CHECK -
                                % Remove profile repetition
                                % (BUG: repeated profiles - SOLVED 20200218)
                                if any(~isnan(platform_irange)) && any(~isnan(JULD_irange))
                                    platform_repeat = ...
                                        ismember(platform(nn),platform_irange);
                                    time_repeat = ...
                                        ismember(JULD(nn),JULD_irange);
                                    
                                else
                                    platform_repeat = 0;
                                    time_repeat = 0;
                                    
                                end
                                
                                rep_prof = sum([platform_repeat time_repeat]);
                                
                                
                                DIST = Distance(lon(nn),lat(nn),lon_sss,lat_sss);
                                irange = find(DIST <= r2);
                                
                                Ln = length(irange);
                                
                                
                                if Ln >= 1 && ismember(rep_prof,[0,1])
                                    ncount = ncount+1;
                                    
                                    % to reset "repetition profile check"
                                    rep_prof = -1;
                                    
                                    %% Satellite data at float location
                                    sss_irange (1:Ln,ncount) = sss(irange);
                                    
                                    sss_error_irange (1:Ln,ncount) = sss_error(irange);
                                    
                                    sss_irange_mn (ncount) = nanmean(sss(irange));
                                    
                                    sss_irange_md (ncount) = nanmedian(sss(irange));
                                    
                                    sss_error_irange_mn (ncount) = nanmean(sss_error(irange));
                                    
                                    sss_error_irange_md (ncount) = nanmedian(sss_error(irange));
                                    
                                    lon_sss_irange (1:Ln,ncount) = lon_sss(irange);
                                    lat_sss_irange (1:Ln,ncount) = lat_sss(irange);
                                    
                                    
                                    %% [2.3] Argo Vertical profile interpolation
                                    [a,b] = size(PSAL);
                                    
                                    xROW = 1:a;
                                    
                                    %nprof_irange (1,ncount) = nprof;
                                    
                                    lon_irange (1,ncount) = lon(nn);
                                    lat_irange (1,ncount) = lat(nn);
                                    
                                    platform_irange (1,ncount) = platform (nn);
                                    ID_irange (1,ncount) = ID(nn);
                                    JULD_irange (ncount) = JULD(nn);
                                    
                                    PRES_irange (xROW,ncount) = PRES(:,nn);
                                    SALT_irange (xROW,ncount) = PSAL(:,nn);
                                    TEMP_irange (xROW,ncount) = TEMP(:,nn);
                                    
                                    clear a b
                                    
                                    %% [2.3.1] Interpolate Z-direction (narrow depth levels)
                                    % Interpolation specs: Grid dimensions
                                    pres1 = 0;              % min pres bin
                                    pres2 = max(PRES);      % max depth (m)
                                    z_grid = 0.5;           % interpolant distance (m) between depth levels
                                    
                                    TEMP1 = TEMP(:,nn);
                                    SALT1 = PSAL(:,nn) ;
                                    PRES1 = PRES(:,nn) ;
                                    
                                    % interpolation function
                                    [TEMP_beta_intp,SALT_beta_intp,PRES_beta_intp] =...
                                        zinterp1_ARGO(TEMP1,SALT1,PRES1,z_grid);
                                    
                                    [a,b] = size(TEMP_beta_intp);
                                    
                                    TEMP_intp(1:a,ncount) = TEMP_beta_intp;
                                    SALT_intp(1:a,ncount) = SALT_beta_intp;
                                    PRES_intp(1:a,ncount) = PRES_beta_intp;
                                    
                                    clear a b *_beta*
                                    
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
            end
            
            
            % compute Matlab time numbe (from Julian to matlab number)
            time_irange = datenum(datetime(JULD_irange,'convertfrom','juliandate'));
            
            %% Remove empty floats (*optimize space)
            % REMEMBER save_vars_out:
            % save_vars_out = [argo_vars_out colocated_vars_out...
            %                stats_X1_vars_out stats_X2_vars_out];
            
            [a,b] = size(platform_irange);
            
            % index empty floats
            
            ind = find(all(isnan(SALT_intp),1));
            
            if ~isempty (ind)
                
                % [BUG] Issue relating to the removing nan (columns/rows) --
                % fixed (badly) with a dirty for-loop and if-statements
                
                for nn = 1:length(save_vars_out)
                    
                    eval(['var_beta = ' save_vars_out{nn} ';']);
                    
                    % To remove platform profiles (NOT depth levels)
                    [a2,b2] = size(var_beta);
                    
                    if b2 == b
                        eval([save_vars_out{nn} '(:,ind) = [];']);
                        
                    elseif a2 == b
                        eval([save_vars_out{nn} '(ind,:) = [];']);
                    end
                    
                end
                
            end;
            clear *beta*
            
            
            %% save fn_out (with Argo-2-satellite colocations)
            
            disp({'New file with Argo-2-satellite colocations '; fn_out})
            save(fn_out,save_vars_out{:})
            
            eval(['clear ' save_vars_out{:}]);
            
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
end


% Setup function output
TT = struct;

for nn = 1:length(save_vars_out)
    this_param = save_vars_out{nn};
    
    eval(['TT.' this_param '= ' this_param ';']);
    
end


vars_out = TT;
clear TT



% END-OF-SCRIPT





