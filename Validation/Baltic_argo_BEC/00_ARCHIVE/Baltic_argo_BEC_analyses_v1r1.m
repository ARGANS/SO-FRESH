% Baltic_argo_BEC_plot
%
% Description
% Analyses of Argo floats as provided by BEC.
% Input: matlab files created from the original netCDF files (use
% ar_rd_BEC.m)
%
% Baltic Study [2011-2013]:
% *[1] Histograms number of floats in region/year
% *[2] Histograms number of floats in region/season
% *[3] Winter Vs Summer profiles (T and S)
%
% Data source:
% BEC products for validation
%
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
% ========================================================================

clc; clear;

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/argo/argo_mat/']);

iyear = 2013; %2011:2012;
imonth = 6:12;

a = length(iyear);
b = length(imonth);

% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file

ibasin = 9; % Basin number: 7 (Arctic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

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
        
        for dd = 1:length(idays)
            iDAY = idays(dd);
            itime_start = datenum(iYEAR,iMONTH,iDAY);
            itime_end = itime_start + ndays;
            
            % load matlab file [fn]
            fn = (['argo_' ...
                datestr(itime_start,'yyyymmdd') '_'...
                datestr(itime_end,'yyyymmdd')]);
            
            fn_in = ([folder_in fn '.mat']);
            
            % ===============
            % Load ARGO float
            if exist(fn_in,'file') == 2
                load (fn_in)
                
                % ===============
                % Load BEC-SSS (neares SSS file to the ARGO file
                if exist(fn_in,'file') == 2
                    load (fn_in)
                    
                elseif exist(fn_in,'file') == 0
                    warning (['file: ' fn_in ' needs to be downloaded'])
                end
                
                % load Satellite SSS
                
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
                
                lon_sss1 = TT.lon;
                lat_sss1 = TT.lat;
                
                % 3/ Baltic-Global product (v001)
                data_type = 6;
                [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,data_type);
                
                sss_gl = TT.sss;
                
                lon_gl = TT.lon;
                lat_gl = TT.lat;
                
                % homogenize grid (global product to baltic+)
                [sss_gl2] = griddata(lon_gl,lat_gl,sss_gl,lon_sss1,lat_sss1);
                
                ind = lon_sss1 >= xmin & lon_sss1 <= xmax |...
                    lat_sss1 >= ymin & lat_sss1 <= ymax ;
                
                sss_gl2 (ind == 0) = NaN;
                
                sss_gl = sss_gl2; clear sss_gl2; % clear work space
                
                
                
                % =================
                % Select SSS data at Float location (within r distance in km)
                
                r = 25;
                % Cases with more than one float (profile) in one file
                
                nprof1 = length(lon(:,1));
                nprof2 = length(lat(:,1));
                
                if nprof1 == nprof2
                    nprof = nprof1;
                    
                else
                    error('Number of Argo longitudes and latitudes must be equal')
                    
                end
                
                
                for nn = length(nprof)
                    DIST = Distance(lon(nn),lat(nn),lon_sss1,lat_sss1);
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
                        
                        % ===============
                        % 1/ boxplot: show to use median of SSS to be compared againts
                        % Argo
                        
                        % string product name
                        PRODS = {'SSS_{NS}','SSS_{NM}','SSS_{gl}'};
                        
                        figure
                        boxplot...
                            ([sss_ns_irange,sss_nm_irange,sss_gl_irange],...
                            'labels',PRODS);
                        
                        ylim([0 10]);
                                                
                        title(['Satellite SSS around Argo float (ID:'...
                            num2str(ID(nn)) ')' ]);
                        
                        
                        
                        
                        % ==============
                        % 2/ Argo profile (T, S)
                        figure
                        
                        
                        
                        
                        
                    end
                end
                
                
            elseif exist(fn_in,'file') == 0
                
                
                % write a log_file, recording the Argo (.mat) files
                file_log = [folder_data '/' 'ARGO_MISSINGt.txt'];
                
                yourMsg = (['file: ' fn_in ' needs to be downloaded']);
                fid = fopen(fullfile(tempdir, file_log), 'a');
                if fid == -1
                    clc
                    warning(yourMsg);
                end
                fprintf(yourMsg);
                %fclose(fid);
            end
        end
    end
end



%% *[1] Plot Argo geolocation
figure (1); clf; hold on

map_projection = 'merc';

lat_min = ymin;
lat_max = ymax;
lon_min = xmin;
lon_max = xmax;

lon_step = 10;
lat_step = 5;

figure; hold on

fillmap_super(map_projection,...
    lon_min,lon_max,lat_min,lat_max,lon_step,lat_step);
hold on


title({['Argo ' basin_str ' (' num2str(lat_min) '-' ...
    num2str(lat_max) '\circ N) ' ]; ...
    [datestr(itime_start,'yyyymmdd') '-'...
    datestr(itime_end,'yyyymmdd')]});
hold on

plotm(lat,lon,'ko','MarkerFaceColor','r')


%% *[2] Plot Argo profiles

%% Interpolate Z-direction (narrow depth levels)
% Interpolation specs: Grid dimensions
pres1 = 0;              % min pres bin
pres2 = max(PRES);      % max depth (m)
z_grid = 0.5;           % interpolant distance (m) between depth levels

[t_inp,s_intp,p_intp] = zinterp1_ARGO(TEMP,SALT,PRES,z_grid);

%% Example
%% 1/ plot raw Vs interp T and S of a given Profile (nprof)
% nprof = find(platform == 6901150); % Location North Atlantic[20161101]
% nprof = find(platform == 4901821); % Location: East Pacific [20161101];

nprof = 1;

if ~isempty(nprof)
    TEMP1 = TEMP (:,nprof);
    SALT1 = PSAL (:,nprof);
    PRES1 = PRES (:,nprof);
    
    TEMP2 = t_intp(:,nprof);
    SALT2 = s_intp(:,nprof);
    PRES2 = p_intp(:,nprof);
    
    plot_argo_BEC(TEMP1,SALT1,PRES1,TEMP2,SALT2,PRES2,platform(nprof))
    
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


