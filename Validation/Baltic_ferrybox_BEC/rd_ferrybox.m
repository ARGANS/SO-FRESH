function [data_out] = rd_ferrybox (fn,varargin)
%
% Syntax: [data_out] = rd_ferrybox (iyear,imonth,iday,ibasin,grid_output)
%
% Output
% data_out = structure with the following content
% data_out = {...
%     'lat','lon',...
%     'SALT','SALT_QC',...
%     'TEMP','TEMP_QC',...
%     'depth','depth_QC',...
%     'date_time', 'date_time_QC'};
%
%
% Description
% Read ferrybox data as provided by BEC for the validation of the Baltic+
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


% path_root = ('/Volumes/Rogue/Data/');
% folder_data = ([path_root ...
%     'SSS/Baltic/BEC/Validation/indata/ferrybox/']);

fn_in = [fn];

% =============
% save vars from netCDF file into the matFILE
param_load = {...
    'LATITUDE','LONGITUDE',...
    'PSAL','PSAL_QC',...
    'TEMP','TEMP_QC',...
    'DEPH','DEPH_QC',...
    'TIME','TIME_QC'};

% variable output names
param_out = {...
    'lat','lon',...
    'SALT','SALT_QC',...
    'TEMP','TEMP_QC',...
    'depth','depth_QC',...
    'date_time', 'date_time_QC'};

param_out = [param_out];

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
    
    
end

time_ref = datenum(1950,01,01,0,0,0);

time_number = nan(size(date_time)); % time matlab number

for nn = 1:length(date_time)
    time_number(nn) = addtodate(time_ref,date_time(nn),'day');
end

% sort obs in chronological order
[time_number,ind] = sort(time_number);

Ln = length(ind); % number of profiles/data

%% [2] Function output (structure)
data_out = struct;



for nn = 1:length(param_out)
    
    eval(['this_param = ' param_out{nn} ';']);
    
    [a,b] = size(this_param);
    
    if  a == Ln  
        
        eval(['data_out. ' param_out{nn} '= ' param_out{nn} '(ind);' ]);
        
    elseif b == Ln
        eval(['data_out. ' param_out{nn} '= ' param_out{nn} '(:,ind);' ]);
        
    end; clear this_param
    
end; clear nn

data_out.time_number = time_number;
    
    
    
end

