function [Output] = ar_rd_BEC_Baltic_FUN(itime,ibasin)
%
% Syntax: [Output] = ar_rd_BEC_Baltic_FUN(iday,imonth,iyear,ibasin)
%
% Input
% itime: matlab time number
% ibasin: basin number
%
% ****************************** ar_nc2mat  ************************************
% Description
% Converted to function: ar_rd_BEC_Baltic.m
% Read NetCDF (.nc) provided by BEC and save them in MATLAB files (.mat).
%
%
%
% --------------------------------------------------------------------------
%
%          Date                     Author                  project
%      10/10/2019               R.E.Jaume-Catany            Baltic+
%
% ----------------------- UPDATES -----------------------------------------
%
% =========================================================================

% Global vars
% ============
% % define the BEC Argo path where data is stored after modifications
path_root = ('/Volumes/Rogue/Data/');
Output = []; % pre-locate output
% ============

% set observing date
[iyear,imonth,iday] = datevec(itime);
iYEARS = iyear;
iMONTHS = imonth;%;1:12;

% ============
% filename Argo BEC structure: argo_20110101_20110110.nc
ndays = 9; % number of days contained in each BEC Argo file

[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);


folder_data = ([path_root ...
    'SSS/' basin_str '/BEC/Validation/indata/argo/']);

% =============
% save vars from netCDF file into the matFILE
param_load = {'CYCLE_NUMBER','DATA_MODE','DAY','DENS','FLAG','ID',...
    'JULD','LATITUDE','LONGITUDE','MONTH','PLATFORM_NUMBER','PRES',...
    'PSAL','RHO','RHO_P0','SANOM','SCLIM','SCLIM','SEC','SU',...
    'TANOM','TCLIM','TEMP','TSU','TU','YEAR','ZSU','ZTU'};

% keep load and save param
param_save = {'CYCLE_NUMBER','DATA_MODE','DAY','DENS','FLAG','ID',...
    'JULD','lat','lon','MONTH','platform','PRES',...
    'PSAL','RHO','RHO_P0','SANOM','SCLIM','SCLIM','SEC','SU',...
    'TANOM','TCLIM','TEMP','TSU','TU','YEAR','ZSU','ZTU'};
% ======
% NetCDF null constants (for more info see ncdisp(fn)) -- Set by BEC
missing_val = -999;
FILL = missing_val;

% =============
% Loop through time to fill-in data
% Define Year and Month
% for yy = 1:length(iYEARS)
%     iyear = iYEARS (yy);
%     for mm = 1: length(iMONTHS)
%         imonth = iMONTHS(mm);

% ===========
% Define folder_in and folder_out
folder_in = ([folder_data 'netcdf/']);
folder_out = ([folder_data 'argo_mat/' ...
    num2str(iyear) '/' sprintf('%02.0f',imonth) '/']);
foldercheck_raf(folder_out);

% ------------------------------------------------------------------------
% DAYS IN A MONTH
% In the next lines I am trying to make that script recognise the length of
% the different folders. This is because each month has different number of
% days. But I can not do it and I just do it manually because time is
% running! (27/04/12)
%         num_days = calendar(iyear,imonth)';
%         num_days(num_days == 0) = [];
%         num_days = length(num_days);

%         for dd = 1 : num_days
%             iday = dd;
date_num1 = datenum(iyear,imonth,iday);
date_num2 = date_num1 + ndays;

filename = (['argo_' ...
    datestr(date_num1,'yyyymmdd') '_'...
    datestr(date_num2,'yyyymmdd')]);

fn_in = ([folder_in filename '.nc']);
fn_out = ([folder_out filename '.mat']);

%                 if exist(fn_in,'file') == 0
%                     error (['file: ' fn_in ' needs to be downloaded'])
%                 end

% =========
% if there is not a matFIle then create it
if exist(fn_in,'file') == 2 && exist(fn_out,'file') ~= 2
    
    % ===========
    % Pre-locate variables (param_save)
    for n = 1 : length(param_save)
        eval([param_save{n} ' = ' num2str(FILL) ';']);
    end
    
    
    % ============
    fill_up = missing_val;  % fill up the file.mat [tmp var]
    
    eval (['save ' fn_out ' fill_up']);
    disp(['I am creating ' fn_out])
    
    param = 'LONGITUDE';
    lon = double(ncread(fn_in,param));
    param = 'LATITUDE';
    lat = double(ncread(fn_in,param));
    
    % Select data within the basin and load the ncFile selection
    ilon = find(lon >= xmin & lon <= xmax);
    ilat = find(lat >= ymin & lat <= ymax);
    
    % Bug#01: nc-files appear to be empty (No data in Baltic)
    % FIX#01: Check lat/lo is not empty.
    if ~isempty(ilon) && ~isempty(ilat)
        % start and counting indices
        start_lon = ilon(1);
        start_lat = ilat(1);
        count_lon = ilon(end)-ilon(1)+1;
        count_lat = ilat(end)-ilat(1)+1;
        
        % load each varriable in the netCDF file
        for x = 1:length(param_load)
            param = param_load{x};
            
            if strcmpi(param,'PLATFORM_NUMBER')
                TT = ncread(fn_in,param);
                TT = str2num(TT');
                eval(['platform = TT;']);
            elseif strcmpi(param,'LONGITUDE')
                TT = double(ncread(fn_in,param));
                ind = TT == missing_val | TT == FILL;
                TT (ind) = NaN;
                eval(['lon = TT;']);
            elseif strcmpi(param,'LATITUDE')
                TT = double(ncread(fn_in,param));
                ind = TT == missing_val | TT == FILL;
                TT (ind) = NaN;
                eval(['lat = TT;']);
            else
                TT = double(ncread(fn_in,param));
                ind = TT == missing_val | TT == FILL;
                TT (ind) = NaN;
                eval([param ' = TT;']);
            end
        end
    end
    % ========
    % ARGO GREY FLOAT
    % Load the grey list (obtained from Argo website.
    % Downloaded Oct 2019: last modified <???>)
    % source: http://www.argo.ucsd.edu/Argo_date_guide.html
    %                     fid = fopen([folder_data 'ar_greylist_20191031.txt']);
    %
    %                     ar_grey = textscan(fid, '%s %s %s %s %s %s %s',...
    %                         'delimiter', ',', 'MultipleDelimsAsOne', 0,'headerlines',1);
    %
    %                     % keep track of the faulty floats
    %                     platform_grey = ar_grey{:,1};   % platform ID
    %                     platform_grey = cellfun(@str2num,platform_grey);
    %                     time_start_grey = ar_grey{:,3}; % time starting being grey
    %                     time_end_grey = ar_grey{:,4};   % time ending being grey
    %
    %                     % index BAD plaforms and remove them
    %                     iBAD = ismember(platform(:),platform_grey(:));
    %
    %                     platform = platform(iBAD);
    
    save(fn_out,param_save{:});
    
    clear (param_save{:});
    
end
disp('work in progresss')


disp(['Argo floats in ' filename '.nc' ...
    ' ' basin_str ' for month: ' num2str(imonth)...
    ' of year '  num2str(iyear)])



% ===========================
% Set the function output

if exist(fn_out,'file') == 2
    Output = load(fn_out);
    
end


