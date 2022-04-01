% Syntax: Baltic_seadatanet_BEC_analyses.m (script)
%
% Description
% Analyses of SeaDataNet data as provided by BEC.
%
%
% 
%
% Input: matlab files created from the original .csv files
%
% Baltic Study [2011-2013]:
% *[1] 
% *[2] 
% *[3]
%
% Data source:
% BEC products for validation
%
%
% current version: v1r0 (2020/01/20)
%
% History
% -
% Version: v1r0
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================

clc; clear;
close all

% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance from reference point (in Km)

ibasin = 9; % Baltic Sea [9];
[xmin,xmax,ymin,ymax,basin_str] = map_lim_raf (ibasin);

fg_save = 1; % flag save figures [1]; or not [0];
fg_format = 'png';

path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/SeaDataNet/']);


% Make a log_file to record status of each ARGO-BEC file [2020/01/21]
folder_log = '/Volumes/Rogue/Data/SSS/Baltic/BEC/Validation/indata/';
fn_log = [folder_log 'SDN_MISSING_20200121.txt'];

folder_figs = ['/Volumes/Rogue/scratch/Validation/SeaDataNet/'];

filename = [folder_data 'data_from_SDN_2015-09_TS_BalticSea_QC_done_v2_filtered.nc'];

iyear = 2013;
imonth = 1:12;


%% [1] Colocation Argo-2-Satellite
run Baltic_SDN_BEC_colocations.m % (snippet)



%% [1.1] Compute Argo at 10 m and the Argo average surface-2-10m
nprof = length(time_irange);

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
    
    while isnan(SALT10_point_beta) && ind_point+1 < zlevels
        ind_point = ind_point+1;
        SALT10_point_beta = SALT_intp(ind_point,nn);
    end
    
    SALT10_point(1,nn) = SALT10_point_beta;
    
end; clear nn *beta*



% function to read SDN dataset
[TT] = rd_SDN (filename);

lon = TT.lon;
lat = TT.lat;

time_number = TT.time_number;
time_str1 = datestr(time_number(1),'yyyymmdd');
time_str2 = datestr(time_number(end),'yyyymmdd');


SALT = TT.SSS_SDN;
TEMP = TT.SST_SDN;
PRES = TT.depth_SDN;

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

%% [1] Plot SDN stations

h1 = plot_map_seadatanet_BEC(lon,lat);

title({...
    ['SeaDataNet ' basin_str]; ...
    [time_str1 '-' time_str2] });

hold on

% Save figure  - output -
fg_name = ['SDN_MAP_' time_str1 '_' time_str2 '.' fg_format];

% Save figure  - output -
folder_this = [folder_figs 'MAPS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];

% check fn existence
fg_exist = exist(fg_name,'file');
if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end



%% [2] Plot Plot Profiles

figure

platform = -999;

nn = 3;
T = TEMP;
S = SALT;
P = PRES;

this_time = time_number(nn);

[h11, h22] = colorplot_argo_platform_profiles(platform,this_time,T,S,P);








