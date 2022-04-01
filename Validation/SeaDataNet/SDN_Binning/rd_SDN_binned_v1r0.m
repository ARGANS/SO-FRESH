function [data_out] = rd_SDN_binned_v1r0 (itime,varargin)
%
% Syntax: [data_out] = rd_SDN_binned_v1r0 (itime)
%
% Description
% Load SeaDataNet data, which has been generated using Baltic_SDN_BEC_Binning_v1r0.
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================

% clc;
close all

REF_str = 'SDN';
iyear = 2012;
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

flag_on = 1;
flag_str = 'FLAG#49';

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

    
    
    
end

