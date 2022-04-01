% Plot SSS Baltic
%
% Description
% Compute monthly climatologies from the daily dataset 2010-2018
%
% Version: v1r0
%
% Date
% 07/06/2019
%
% rcatany
%
%
% =======================================================================
%DATA_PATH = '/Users/rejc1e11/Documents/'; % # Use this path working local
DATA_PATH = '/Volumes/Rogue/Data/';            % # Use this path working hard drive
DATA_PRODUCT = 'global';                    % # BEC product type: [1] Global; [2] Arctic (avail. Sep 2019); [3] Baltic (avail Sep 2019)

% Srcipt settings
ibasin = 9; % Basin number (Baltic, ibasin 9)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);


iYEAR = 2011:2016;      % years range
iMONTH = 1:12;          % month range
iday = 1;

nYEAR = length(iYEAR);          % number of years
nMONTH = length(iMONTH);        % number of months

% min/max SSS
Smin = 0;
Smax = 35;

% ==========================

% Matrix pre-location
[A] = rd_smos_L4_BEC_v1r1(datenum(iYEAR(1),iMONTH(1),1),ibasin); % read any sss to get dimensions (a, b)

sss = A.sss;
lon = A.lon;
lat = A.lat;


[a,b] = size(sss);
[c] = 31; % maximum number of days in a month
[d] = nMONTH;
[e] = nYEAR;

sss_month = NaN(a,b,c);

% ==============
% Construct monthly averages (save files)

for yy = 1:nYEAR
    for mm = 1:nMONTH
        iyear = iYEAR(yy);
        imonth = iMONTH (mm);
        DAYS = calendar(iyear,imonth); DAYS = sort(DAYS(DAYS~=0));
        nDAYS = length(DAYS);
        
        %
        if ibasin == 7
            folder_out = [DATA_PATH 'SSS/Arctic/BEC/monthly/' num2str(iyear) '/'];
            
        elseif ibasin ~= 7
            folder_out = [DATA_PATH 'SSS/SMOS/BEC/data/global/monthly/' num2str(iyear) '/'];
        end
        
        foldercheck_raf (folder_out);
        
        fn_out = [folder_out basin_str '_MONTHLY_' sprintf('%02.0f',iyear)...
            '_' sprintf('%02.0f',imonth) '.mat'];
        
        if exist(fn_out,'file') ~= 2
            
            % # Inform about what month and year you are looping on.
            disp(['Get month averagles for month ' ...
                datestr(datenum(iyear,imonth,1),'yyyy-mm')]);
            
            for dd = 1:nDAYS
                iday = DAYS(dd);
                
                [A] = rd_smos_L4_BEC_v1r1(datenum(iyear,imonth,iday),ibasin);
                sss_month (:,:,dd) = A.sss;
                clc
            end % # end daily loop
            
            sss_mean = nanmean(sss_month,3);
            sss_median = nanmean(sss_month,3);
            
            % # Save fn_out with monthly averages (mean, median)
            save (fn_out, 'sss_mean', 'sss_median','lon','lat');
        end
        
    end % # end monthly loop [yy]
end % # end yearly loop [mm]

clear sss_mean sss_median sss_month yy mm fn_out

% ============
%% Construct Climatology (2011-2017)

[a,b] = size(sss);
[c] = nMONTH;
[d] = nYEAR;

SSS_month_mean = NaN(a,b,c,d);
SSS_month_median = NaN(a,b,c,d);

% 1) Get each month in every year (EMEY)
for mm = 1:nMONTH
    for yy = 1:nYEAR
        iyear = iYEAR(yy);
        imonth = iMONTH (mm);
        
        folder_in = [DATA_PATH 'SSS/SMOS/BEC/data/' DATA_PRODUCT '/monthly/'...
            num2str(iyear) '/' ];
        
        fn_in = [folder_in basin_str '_MONTHLY_' num2str(iyear) '_'...
            sprintf('%02.0f',imonth) '.mat'];
        
        load(fn_in,'sss_mean','sss_median','lon','lat')
        
        SSS_month_mean(:,:,mm,yy) = sss_mean;
        SSS_month_median(:,:,mm,yy) = sss_median;
        
        clear sss_mean sss_median
        
    end % # get EMEY
end % # end monthly climatology [mm]


% # 2) Average each month with the year range (2011-2017)

folder_out = [DATA_PATH 'SSS/SMOS/BEC/data/' DATA_PRODUCT '/Climatology/'];
foldercheck_raf (folder_out);

for mm = 1:nMONTH
    imonth = iMONTH (mm); 
    
    fn_out = [folder_out basin_str '_CLIMAT_' num2str(iYEAR(1)) num2str(iYEAR(end)) ...
        '_' sprintf('%02.0f',imonth) '.mat'];
    
    if exist(fn_out,'file') ~= 2
    
    SSS_climat_mean = squeeze(SSS_month_mean(:,:,mm,:));
    SSS_climat_median = squeeze(SSS_month_mean(:,:,mm,:));
    
    SSS_climat_mean (:,:,mm) = mean(SSS_climat_mean,3);
    SSS_climat_median (:,:,mm) = median(SSS_climat_median,3);
    
    save(fn_out,'SSS_climat_mean','SSS_climat_median', 'lon','lat');
    
    end
    
end









