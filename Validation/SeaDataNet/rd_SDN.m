function [data_out] = rd_SDN (fn,varargin)
%
% Syntax: [data_out] = rd_SDN (iyear,imonth,iday,ibasin,grid_output)
%
% Description
% Read SeaDataNet data as provided by BEC for the validation of the Baltic+
% product.
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================

% switch off warning (there is one about addtodate function)
warning('off','all')

if isempty (varargin)
    griddata_on = 0;
else
    griddata_on = 1;
end


% Get SSS-satellite data within a given distance from each Argo float
r = 25; % Radius distance from reference point (in Km)


path_root = ('/Volumes/Rogue/Data/');
folder_data = ([path_root ...
    'SSS/Baltic/BEC/Validation/indata/SeaDataNet/']);

% Sugested filename convention
% fn = ['SDN_YYYYMMDD_YYYYMMDD_TS_BalticSea_QC_v2_filtered'];


fn_in = [fn];

% =============
% save vars from netCDF file into the matFILE
param_load = {'latitude','longitude',...
    'var1','var1_qc',...
    'var2','var2_qc',...
    'var3','var3_qc',...
    'date_time'};

metavar_load = {'metavar3','metavar5','metavar12'};

% variable output names
param_out = {'lat','lon',...
    'depth_SDN','depth_QC',...
    'SST_SDN','SST_QC',...
    'SSS_SDN','SSS_QC',...
    'date_time'};

param_out = [param_out metavar_load];

% NetCDF null constants (for more info see ncdisp(fn)) -- Set by BEC
missing_val = -999;
FILL = missing_val;



%% [1] read each variable in the NetCDF file
if exist(fn_in,'file') == 2
    
    % [1] Load parameters (numeric values)
    for n = 1:length(param_load)
        param = param_load{n};
        TT = squeeze(double(ncread(fn_in,param)));
        ind = TT == missing_val | TT == FILL;
        TT (ind) = NaN;
        eval([param_out{n} ' = TT;']);
        
        clear TT;
    end; clear n
    
    
    % [2] Load metadata (string with instrument informaiton)
    for nn = 1:length(metavar_load)
        param = metavar_load{nn};
        TT = squeeze(char(ncread(fn_in,param)));
        ind = TT == missing_val | TT == FILL;
        TT (ind) = NaN;
        eval([metavar_load{nn} ' = TT;']);    
              
        
        clear TT;
    end; clear nn
end

time_ref = datenum(2011,01,01,0,0,0);

time_number = nan(size(date_time)); % time matlab number

for nn = 1:length(date_time)
    time_number(nn) = addtodate(time_ref,date_time(nn),'day');
end

% sort obs in chronological order
[time_number,ind] = sort(time_number);

Ln = length(ind); % number of profiles/data

%% [2] Function output (structure)
data_out = struct;

data_out.time_number = time_number;

for nn = 1:length(param_out)
    
    eval(['this_param = ' param_out{nn} ';']);
    
    [a,b] = size(this_param);
    
    if  a == Ln  
        eval(['data_out. ' param_out{nn} '= ' param_out{nn} '(ind);' ]);
        
    elseif b == Ln
        eval(['data_out. ' param_out{nn} '= ' param_out{nn} '(:,ind);' ]);
        
    end; clear this_param
    
end; clear nn


    
    
    
end

