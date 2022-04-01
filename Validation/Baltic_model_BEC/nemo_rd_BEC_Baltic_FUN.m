function [data_model] = nemo_rd_BEC_Baltic_FUN (iyear,imonth,iday,ibasin,varargin)
%
% Syntax: [data_model] = nemo_rd_BEC_Baltic_FUN (iyear,imonth,iday,ibasin,grid_output)
%
% Description
% Read 3D model for Baltic+ Validation . The model is:
% "CMEMS V4 Reanalysis: NEMO model 3D fields (monthly means)"
%
%  >> For the validation purposes of the Baltic+ there are two periods worth of
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
% Input
% [1] iyear : from 2011 to 2016
% [2] imonth 
% [3] iday
% [4] ibasin (9: for Baltic limits)
% [5] Griddata [1] output [optional] to SMOS-BEC (Baltic+) grid, or not [0]
%
% Note: Set iday = 0 and imonth = 0 To get all data in each file.
%
% Output 
% data structure with the following content:
%
% tt.SSS_model = squeeze(SSS_model(:,:,ind));
% tt.SST_model = squeeze(SST_model(:,:,ind));
% tt.MLD_model = squeeze(MLD_model(:,:,ind));
% tt.lon = lon;
% tt.lat = lat;
% tt.TIME_model = TIME_model;
% 
% =========================================================================
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

% remove matlab warning (setting off the conversion matlab number)
warning('off','all') 

if isempty (varargin)
    griddata_on = 0;
else
    griddata_on = 1;
end

itime = datenum(iyear,imonth,iday);

% to choose 2011 or 2014
if iyear >= 2011 && iyear <= 2013
    year1 = 2011;
elseif iyear >= 2014 && iyear <= 2016
    year1 = 2014;
else
    error(['Input: ' datestr(itime,'yyyymmdd')...
        ' out of temporal limits: 2011-2013, 2014-2016'])
end

year2 = year1+2;

year1_str = num2str(year1);
year2_str = num2str(year2);


[~,~,~,~,basin_str] = map_lim_raf (ibasin);

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


% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'NEMO_MISSING_20200121.txt'];

fn = ['dataset-reanalysis-nemo-monthlymeans_' year1_str '_' year2_str];
folder_in = [folder_data];
fn_in = [folder_in fn '.nc'];

fn_exist = exist(fn_in,'file');

if fn_exist ~= 2
    error(['To Check your temporal input, '...
        'Baltic-NEMO data is only avalable within 2011 to 2016'])
end


% =============
% save vars from netCDF file into the matFILE
param_load = {'latitude','longitude','mlotst','so','thetao','time','depth'};

% NetCDF null constants (for more info see ncdisp(fn)) -- Set by BEC
missing_val = -999;
FILL = missing_val;

%% [1] read each variable in the NetCDF file
if exist(fn_in,'file') == 2
    for x = 1:length(param_load)
        param = param_load{x};
        if strcmpi(param,'so')
            TT = squeeze(double(ncread(fn_in,param)));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['SSS_model = TT;']);
        elseif strcmpi(param,'thetao')
            TT = squeeze(double(ncread(fn_in,param)));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['SST_model = TT;']);
        elseif strcmpi(param,'LONGITUDE')
            TT = squeeze(double(ncread(fn_in,param)));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['lon = TT;']);
        elseif strcmpi(param,'LATITUDE')
            TT = squeeze(double(ncread(fn_in,param)));
            ind = TT == missing_val | TT == FILL;
            TT (ind) = NaN;
            eval(['lat = TT;']);
        elseif strcmpi(param,'mlotst')
            TT = squeeze(double(ncread(fn_in,param)));
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




% =============================
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

% If especific date is not selected, keep all data output
if iday ~= 0 && imonth ~= 0
    
    ind = find(itime >= time_out,1);
    TIME_model = time_out(ind);
    
else
    TIME_model = time_out;
    ind = 1:length(TIME_model);
end


% =============================
% Convert NEMO-grid to SMOS-BEC grid


if griddata_on == 1
    
    folder_out = folder_in;
    fn_out = [folder_out fn '_GRIDDED_025.mat'];
    
    fn_exist = exist(fn_out,'file');
    
    save_vars = {'SSS_model','SST_model','MLD_model','TIME_model','lon','lat'};
    
    if fn_exist ~= 2
        clc
        warning (['I am gridding NEMO output to SSS-BEC Baltic+. '...
            'It will take a few secs. '])
        
        disp(['Output file: ' fn_out ])
        
        
        data_type = 4;
        
        itime_smos = datenum(2011,09,01,0,0,0); % selected a radom file to get lon/lat
        
        [A] = rd_smos_BEC_Baltic_FUN(itime_smos,ibasin,data_type);
        lon_smos = A.lon;
        lat_smos = A.lat;
        clear A
        
        
        [a,b] = size(lon_smos);
        c = length(TIME_model);
        
        SSS_model_gridded = NaN(a,b,c);
        SST_model_gridded = NaN(a,b,c);
        MLD_model_gridded = NaN(a,b,c);
        
        
        for nn = 1:length(TIME_model)
            SSS_model_gridded (:,:,nn) = griddata(lon,lat,SSS_model(:,:,nn),lon_smos,lat_smos);
            SST_model_gridded (:,:,nn) = griddata(lon,lat,SST_model(:,:,nn),lon_smos,lat_smos);
            MLD_model_gridded (:,:,nn) = griddata(lon,lat,MLD_model(:,:,nn),lon_smos,lat_smos);
        end
        
        
        lon = lon_smos;
        lat = lat_smos;
        
        SSS_model = SSS_model_gridded;
        SST_model = SST_model_gridded;
        MLD_model = MLD_model_gridded;
        
        clear *_smos *_gridded
        
        save(fn_out,save_vars{:})
        
    else
        load (fn_out)
    end
    
    
end

% =============================
%% [2] Statistics in each Baltic-region
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



%% [2.1] Get NEMO-data within a searching distance at each study region
% Get NEMO-data within a given distance from the center if SREG.
r = 200; % km

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

lon_irange = nan(a*b,Nsreg);
lat_irange = nan(a*b,Nsreg);

for nn = 1:length(lon_sreg)
    DIST = Distance(lon_sreg(nn),lat_sreg(nn),lon,lat);
    
    irange = find(DIST <= r);
    Ln = length(irange);
    
    lon_irange(1:Ln,nn) = lon(irange);
    lat_irange(1:Ln,nn) = lat(irange);
    
    SSS_irange(1:Ln,:,nn) = SSS_model_draft(irange,:);
    SST_irange(1:Ln,:,nn) = SST_model_draft(irange,:);
    MLD_irange(1:Ln,:,nn) = MLD_model_draft(irange,:);
    
end


% =======================================
% Return to original shape NEMO-inrange variables
% [a,b,c] = size(SSS_model);
%
% SSS_irange = reshape(SSS_irange,a,b,c,Nsreg);
% SST_irange = reshape(SST_irange,a,b,c,Nsreg);
% MLD_irange = reshape(MLD_irange,a,b,c,Nsreg);
%
% lon_irange = reshape(lon_irange,a,b,Nsreg);
% lat_irange = reshape(lat_irange,a,b,Nsreg);


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



% ============================
% Keep function structure output
tt = struct;
tt.SSS_model = squeeze(SSS_model(:,:,ind));
tt.SST_model = squeeze(SST_model(:,:,ind));
tt.MLD_model = squeeze(MLD_model(:,:,ind));
tt.lon = lon;
tt.lat = lat;
tt.TIME_model = TIME_model;

% Regional Baltic data and statistics

tt.lat_sreg = lat_sreg;
tt.lon_sreg = lon_sreg;

tt.lat_irange = lat_irange;
tt.lon_irange = lon_irange;

tt.SREG = SREG;
tt.SREG_str = SREG_str;
tt.SREG_str2 = SREG_str2;

tt.SSS_irange = SSS_irange;
tt.SST_irange = SST_irange;
tt.MLD_irange = MLD_irange;

tt.SSS_mean_sreg = SSS_mean_sreg;
tt.SST_mean_sreg = SST_mean_sreg;
tt.MLD_mean_sreg = MLD_mean_sreg;

tt.SSS_median_sreg = SSS_median_sreg;
tt.SST_median_sreg = SST_median_sreg;
tt.MLD_median_sreg = MLD_median_sreg;

tt.SSS_std_sreg = SSS_std_sreg;
tt.SST_std_sreg = SST_std_sreg;
tt.MLD_std_sreg = MLD_std_sreg;

tt.SSS_mean_sreg_ALL = SSS_mean_sreg_ALL;
tt.SST_mean_sreg_ALL = SST_mean_sreg_ALL;
tt.MLD_mean_sreg_ALL = MLD_mean_sreg_ALL;

tt.SSS_median_sreg_ALL = SSS_median_sreg_ALL;
tt.SST_median_sreg_ALL = SST_median_sreg_ALL;
tt.MLD_median_sreg_ALL = MLD_median_sreg_ALL;

tt.SSS_std_sreg_ALL = SSS_std_sreg_ALL;
tt.SST_std_sreg_ALL = SST_std_sreg_ALL;
tt.MLD_std_sreg_ALL = MLD_std_sreg_ALL;



data_model = tt; clear tt ind

% ============================

