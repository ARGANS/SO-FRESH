function Baltic_plot_fun(ibasin,itime,data_type,varargin)
%
% Syntax: Baltic_plot_fun(ibasin,itime,data_type,plot_type,region)
%
% Description
% Plot to see SSS BEC products in the Arctic.
% This function builds up on the script Arctic_plot_v1r1.m
%
% Input (as seen in function rd_smos_L4_BEC_v1r2.m)
% itime: date as in MATLAB number
% ibasin: basin number (for more info type help map_lim_raf)
%
% data_type for Baltic:
%           [1] 9_days Arctic+ product (v3.0)
%           [2] 9_days Arctic BEC product (v2.0)
%           [3] Global BEC product (v001)
%           [4] Baltic BEC NOM product (v1.0)
%           [5] Baltic BEC NS product (v1.0)
%           [6] Global (default)
%           [7] CCI+SSS
%
% BEC_PRODUCTS = [1:6];
% plot_type for Arctic
%           [1] Mercartor
%           [2] Polar
%           [3] Lambert (conical) [plot each side of the Arctic: Atlantic, Pacific
%           sides]
%           [4] All- Plot all projection above
% plot_type for Baltic
%           [1] Mercator
%           [2] Lambert (conical)
% region: 'Baltic'
%
% Version: v1r0
% Author: rcatany (rafaelcatany@gmail.com)
%
% History
% Version | Date| Note
% v1r0 [24/10/2019] : Creation of the function
% =======================================================================
%%

% clc;
% clear

% Example (settings)
% ibasin = 7;     % Basin number (Arctic, ibasin 7)
% data_type: [1-3] Arctic+
% ([1]monthly,[2] 7day or [3] 9 day), Baltic [4,5], product and [6] Global product
% data_type = 3;
%
% iyear = 2014;
% imonth = 2;
% iday = 1;
%
% itime = datenum(iyear,imonth,iday);


DATA_PATH = '/Volumes/Rogue/Data/';            % # Use this path working hard drive


if nargin == 3
    plot_type = 2;  % Mercator projection *set as default
    region = 'ALL'; % make all region * set as default
    
elseif nargin == 4
    plot_type = varargin{1};
    region = 'ALL';
    
elseif nargin == 5
    plot_type = varargin{1};
    region = varargin{2};
    
end

[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

save_fig = 1; % flag to save figure [01]; or not [00].
fg_format = 'png';

folder_out = ['/Volumes/Rogue/scratch/MAPS/' basin_str]; % folder to store figures

% BEC_PRODUCTS flagged from 1 to 6
BEC_PRODUCTS = [1:6];


if data_type == 1
    % ARC_PRODUCT = ibasin == 7 &&  data_type == 1;% Arctic+ product (v3.0)
    data_type_str = ['ARCTIC+ product (v3.0)'];
    data_version = 'ARCTIC-BECv3r0';
    
elseif data_type == 2
    % ARC_PRODUCT = ibasin == 7 && data_type == 2; % Arctic BEC prodcut (v2.0)
    data_type_str = ['ARCTIC product (BEC v2.0)'];
    data_version = 'ARCTIC-BECv2r0';
    
elseif data_type == 4
    data_type_str = ['Baltic Nominal product (v1.0)'];
    data_version = 'BALTIC-BEC-NMv1r0';
    
elseif data_type == 5
    data_type_str = ['Baltic NS product (v1.0)'];
    data_version = 'BALTIC-BEC-NSv1r0';
elseif data_type == 6
    % ARC_PRODUCT = 0;
    data_type_str = ['GLOBAL product (v001)'];
    data_version = 'GLOBAL-BECv001';
    
elseif data_type == 7
    data_type_str = ['GLOBAL product CCI+SSS (v1.8)'];
    data_version = 'CCI+SSSv1.8';    
    
end


% create a storage folder
folder_out = [folder_out '/' data_version];


if save_fig == 1
    foldercheck_raf(folder_out); %! make folder_out
    cd (folder_out)
    !open pwd
end


% ============
% Salinity limits (plot)

if ibasin == 7
    Smin = 31;
    Smax = 37;
    S_error_min = 0;
    S_error_max = 2;
    
elseif ibasin == 9
    Smin = 0;
    Smax = 10;
    S_error_min = 0;
    S_error_max = 2;
end

% ==========================
% load satellite product as BEC_PRODUCTS, or a CCI+SSS product


if ismember(data_type,BEC_PRODUCTS)%
    
    [A] = rd_smos_L4_BEC_v1r2(itime,ibasin,data_type);
    
elseif data_type == 7 % CCI+SSS product
    [A] = rd_sss_cci(itime,ibasin,1);
    
end

lon  = A.lon;
lat = A.lat;
sss = A.sss;



% ==================================


%% load Baltic+ grid -to convert all products to the same grid
GLOBAL_PRODUCT = [6,7];

if ismember(data_type, GLOBAL_PRODUCT)
    
    % store Baltic+ grid (lon/lat)
    folder_grid = ([DATA_PATH ...
        'SSS/' basin_str '/BEC/Baltic_grid/']);
    
    fn_grid = [folder_grid 'Baltic_plusv1.0_grid.mat'];
    
    if exist(fn_grid,'file') ~= 2
        [TT] = rd_smos_L4_BEC_v1r3(itime_start,ibasin,4);
        
        lon_sss = TT.lon;
        lat_sss = TT.lat;
        
        % save Arctic grid
        folder_grid = (['/Volumes/Rogue/Data/SSS/Baltic/BEC/' basin_str '_grid/']);
        foldercheck_raf(folder_grid);
        
        save (fn_grid,'lon_sss','lat_sss');
        
    else
        load (fn_grid)
    end
    
    
    %% Keep the gridded output (i.e. same grid as Baltic+ product)
    sss_alfa = griddata(lon,lat,sss,lon_sss,lat_sss);

    % remove data from outside the Baltic sea (outside Arkona basin)
    lon_sss(lon_sss < 14) = NaN;
    lat_sss(lat_sss > 66.5) = NaN;
    
    
    lon = lon_sss;
    lat = lat_sss;
    sss = sss_alfa;
    
    clear *_sss _alfa
    
end


% flag baltic+ and arctic+
salinity_plus = [1,4,5];

if ismember(data_type, salinity_plus)
    sss_error = A.sss_error;
end


time_out = A.time_out;

clear A; clc;

A = size(lon);
B = size(lat);
C = size(sss);

if isequal(A,B) == 0 % lon/lat in Baltic product is not equal: nxm
    [a,b] = size(sss);
    lat = repmat(lat',[90 1]);
    lon = repmat(lon,[1 58]);
end


% lon_min = -180;
% lon_max = 180;
% lat_min = 40;
% lat_max = 90;

lon_min = xmin;
lon_max = xmax;
lat_min = ymin;
lat_max = ymax;


[IND] = find(lon >= lon_min & lon <= lon_max &...
    lat >= lat_min & lat <= lat_max);

lon2 = lon(IND);
lat2 = lat(IND);
sss2 = sss(IND);

xLON = lon_min+0.25:0.25:lon_max-0.25;
yLAT = lat_min-0.25:0.25:lat_max-0.25;

[xLON2,yLAT2] = meshgrid(xLON,yLAT);

% Keep the gridded output (lon3,lat3,sss3)
sss3 = griddata(lon2,lat2,sss2,xLON2,yLAT2);


% Load and Grid SSS error (only for Salinity+ data
if ismember(data_type, salinity_plus)
    sss3_error = sss_error(IND);
    sss3_error = griddata(lon2,lat2,sss3_error,xLON2,yLAT2);
end

lon3 = xLON2;
lat3 = yLAT2;

% # =======================================
%% 01) Plot the contentent of the Arctic dataset (Mercator projection)

% Salinity limits (plot)
% Smin = 25;
% Smax = 38;
if plot_type == 1
    % Set the map projection specifications
    map_projection = 'merc';
    lon_min = xmin;
    lon_max = xmax;
    lat_min = ymin;
    lat_max = ymax;
    
    lon_step = 10;
    lat_step = 5;
    
    figure; hold on
    
    fillmap_super(map_projection,...
        lon_min,lon_max,lat_min,lat_max,lon_step,lat_step);
    hold on
    pcolorm(lat3,lon3,sss3); shading flat
    
    colormap(jet)
    caxis ([Smin Smax])
    colorbar
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1].*0.5)
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue')
    gridm on
    title({[basin_str ' SSS ' datestr(itime,'yyyymmdd')];[data_type_str]})

    % Save figure  - output -
   
    fg_name = [data_version ...
        '_' 'PROJ#MERC' '_' datestr(itime,'yyyymmdd') '_SSS.' fg_format];
    
    % store figs in this folder
    folder_out_this = [folder_out '/SSS'];    
    foldercheck_raf(folder_out_this); %! make folder_out

    fg_name = [folder_out_this '/' fg_name];
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if save_fig == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
    
    % Make a map of SSS error (only for Salinity+ data
    if ismember(data_type, salinity_plus)
        
        figure; hold on
        
        fillmap_super(map_projection,...
            lon_min,lon_max,lat_min,lat_max,lon_step,lat_step);
        hold on
        pcolorm(lat3,lon3,sss3_error); shading flat
        
        colormap(jet)
        caxis ([S_error_min S_error_max])
        colorbar
        framem on
        axis off
        tightmap
        geoshow('landareas.shp','facecolor',[1 1 1].*0.5)
        rivers = shaperead('worldrivers', 'UseGeoCoords', true);
        geoshow(rivers, 'Color', 'blue')
        gridm on
        title({[basin_str ' SSS_{error} ' ...
            datestr(itime,'yyyymmdd')];[data_type_str]})
        
        % Save figure  - output -
        fg_name = [data_version ...
            '_' 'PROJ#MERC' '_' datestr(itime,'yyyymmdd') '_SSSerror.' fg_format];
        
        % store figs in this folder
        folder_out_this = [folder_out '/SSSerror'];
        foldercheck_raf(folder_out_this); %! make folder_out
        
        fg_name = [folder_out_this '/' fg_name];
        
        fg_exist = exist(fg_name,'file'); % check fn existence
        
        if save_fig == 1 && fg_exist == 0
            save_raf(gcf,fg_name,fg_format); close
        end
         
        
    end
end
% ========================================
%%
%% 03) Polar plot projection (follow matlab map limit properties
%% and axesm documentation)
%
% Salinity limits (plot)
% Smin = 25;
% Smax = 38;
if plot_type == 2
    load coastlines
    
    % [fig 2]
    figure; clf; hold on
    axesm('eqaazim','MapLatLimit',[50 90])
    axis off
    framem on
    gridm on
    mlabel on
    plabel on;
    setm(gca,'MLabelParallel',40)
    
    pcolorm(lat3,lon3,sss3); shading flat
    colormap(jet)
    caxis ([Smin Smax])
    colorbar
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1].*0.5)
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue')
    gridm on
    title({[basin_str ' SSS ' datestr(itime,'yyyymmdd')];[data_type_str]})
    
    % Save figure  - output -
    fg_name = [data_version ...
        '_' 'PROJ#POLAR' '_' datestr(itime,'yyyymmdd') '_SSS.' fg_format];
    
    % store figs in this folder
    folder_out_this = [folder_out '/SSSerror'];
    foldercheck_raf(folder_out_this); %! make folder_out
    
    fg_name = [folder_out_this '/' fg_name];
    
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if save_fig == 1 && fg_exist == 0
        
        save_raf(gcf,fg_name,fg_format); close
        
    end 
    
    
end


%% 03) Lambert plot projection (Conical)
% Selection Arctic regions of interest (as seen in DUM?)
%
% 1* Subpolar North Atlantic (40?N 65?N, 90?W 30?E)
% 2* Polar North Atlantic
% 3* North West Arctic Pacific
% 4* North East Arctic Pacific

% 1* Subpolar North Atlantic (40?N 65?N, 90?W 30?E)

if plot_type == 3 && strcmpi(region,'NA') || strcmpi(region,'ALL')
    
    % [fig 3]
    figure; clf
    lat_min2 = 50;
    lat_max2 = 85;
    lon_min2 = -100;
    lon_max2 = 100;
    
    lat_step = 15;
    lon_step = 20;
    
    Smin2 = 25;
    Smax2 = 38;
    
    map_projection = 'lambert';
    
    fillmap_super(map_projection,lon_min2,lon_max2,...
        lat_min2,lat_max2,lon_step,lat_step)
    
    pcolorm(lat3,lon3,sss3); shading flat
    colormap(jet)
    caxis ([Smin2 Smax2])
    colorbar
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1]*.5)
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue')
    gridm on
    title({[basin_str ' North Atlantic (' num2str(lat_min2) ...
        '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
        [data_type_str]})
    
    % Save figure  - output -
    fg_name = [data_version ...
        '_' 'PROJ#LAMB' '_' datestr(itime,'yyyymmdd') '_SSS.' fg_format];
    
    fg_name = [folder_out '/' fg_name];
    
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if save_fig == 1 && fg_exist == 0
        
        save_raf(gcf,fg_name,fg_format); close
        
    end 
    
    
end

%% 2* Polar North Atlantic (60N 90N, 0 180E)

if plot_type == 3 && strcmpi(region,'NA') || strcmpi(region,'ALL')
    
    % [fig 4]
    figure; clf
    lat_min2 = 60;
    lat_max2 = 85;
    lon_min2 = 0;
    lon_max2 = 180;
    
    lat_step = 15;
    lon_step = 20;
    
    Smin2 = 33;
    Smax2 = 36;
    
    map_projection = 'lambert';
    
    fillmap_super(map_projection,lon_min2,lon_max2,...
        lat_min2,lat_max2,lon_step,lat_step)
    
    
    pcolorm(lat3,lon3,sss3); shading flat
    colormap(jet)
    caxis ([Smin2 Smax2])
    colorbar
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1]*.5)
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue')
    gridm on
    title({[basin_str ' North Atlantic (' num2str(lat_min2) ...
        '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
        [data_type_str]})
    
    % Save figure  - output -
    fg_name = [data_version ...
        '_' 'PROJ#LAMB' '_' datestr(itime,'yyyymmdd') '_SSS.' fg_format];
    
    % store figs in this folder
    folder_out_this = [folder_out '/SSS'];
    foldercheck_raf(folder_out_this); %! make folder_out
    
    fg_name = [folder_out_this '/' fg_name];
    
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if save_fig == 1 && fg_exist == 0
        
        save_raf(gcf,fg_name,fg_format); close
        
    end 
    
    
end

%% 3* North West Arctic Pacific
if plot_type == 3 && strcmpi(region,'PA') || strcmpi(region,'ALL')
    % [fig 5]
    figure; clf
    lat_min2 = 45;
    lat_max2 = 80;
    lon_min2 = 130;
    lon_max2 = -120;
    
    lat_step = 15;
    lon_step = 20;
    
    Smin2 = 30;
    Smax2 = 34;
    
    map_projection = 'lambert';
    
    
    fillmap_super(map_projection,lon_min2,lon_max2,...
        lat_min2,lat_max2,lon_step,lat_step)
    
    setm(gca, 'Origin', [0 -180]) % Set the Pacific in the middle of the plot
    
    pcolorm(lat,lon,sss); shading flat
    colormap(jet)
    caxis ([Smin2 Smax2])
    colorbar
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1]*.5)
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue')
    gridm on
    title({[basin_str ' North Pacific (' num2str(lat_min2) ...
        '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
        [data_type_str]})
    
    % Save figure  - output -
    fg_name = [data_version ...
        '_' 'PROJ#LAMB' '_' datestr(itime,'yyyymmdd') '_SSS.' fg_format];
    
    % store figs in this folder
    folder_out_this = [folder_out '/SSS'];
    foldercheck_raf(folder_out_this); %! make folder_out
    
    fg_name = [folder_out_this '/' fg_name];
    
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if save_fig == 1 && fg_exist == 0
        
        save_raf(gcf,fg_name,fg_format); close
        
    end
    
    
end


%% 03) Lambert plot projection (Conical)
% Selection Baltic regions of interest (as seen in DUM?)
%
% 1* Subpolar North Atlantic (40?N 65?N, 90?W 30?E)
% 2* Polar North Atlantic
% 3* North West Arctic Pacific
% 4* North East Arctic Pacific


if plot_type == 3 && strcmpi(region,'Baltic') ...
        || data_type == 4 && strcmpi(region,'Baltic')...
        || data_type == 5 && strcmpi(region,'Baltic')
    
    % Plot SSS map
    figure; clf
    lat_min2 = ymin;
    lat_max2 = ymax;
    lon_min2 = xmin;
    lon_max2 = xmax;
    
    lat_step = 2.5;
    lon_step = 4;
    
    
    map_projection = 'lambert';
    
    fillmap_super(map_projection,lon_min2,lon_max2,...
        lat_min2,lat_max2,lon_step,lat_step)
    
    pcolorm(lat3,lon3,sss3); shading flat
    colormap(jet)
    caxis ([Smin Smax])
    colorbar
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1]*.5)
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue')
    gridm on
    title({[basin_str ' (' num2str(lat_min2) ...
        '-' num2str(lat_max2) '\circ N) SSS ' datestr(itime,'yyyymmdd')]; ...
        [data_type_str]})
    
    
    % Save figure  - output -
    fg_name = [data_version ...
        '_' 'PROJ#LMB' '_' datestr(itime,'yyyymmdd') '_SSS.' fg_format];
    
    % store figs in this folder
    folder_out_this = [folder_out '/SSS'];
    foldercheck_raf(folder_out_this); %! make folder_out
    
    fg_name = [folder_out_this '/' fg_name];
    
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if save_fig == 1 && fg_exist == 0
        
        save_raf(gcf,fg_name,fg_format); close
        
    end
    
    
    
    % Make a map of SSS error (only for Salinity+ data
    if ismember(data_type, salinity_plus)
        figure; clf
        lat_min2 = ymin;
        lat_max2 = ymax;
        lon_min2 = xmin;
        lon_max2 = xmax;
        
        lat_step = 2.5;
        lon_step = 4;
        
        
        map_projection = 'lambert';
        
        fillmap_super(map_projection,lon_min2,lon_max2,...
            lat_min2,lat_max2,lon_step,lat_step)
        
        pcolorm(lat3,lon3,sss3_error); shading flat
        colormap(jet)
        caxis ([S_error_min S_error_max])
        colorbar
        framem on
        axis off
        tightmap
        geoshow('landareas.shp','facecolor',[1 1 1]*.5)
        rivers = shaperead('worldrivers', 'UseGeoCoords', true);
        geoshow(rivers, 'Color', 'blue')
        gridm on
        title({[basin_str ' (' num2str(lat_min2) ...
            '-' num2str(lat_max2) '\circ N) SSS_{error} ' datestr(itime,'yyyymmdd')]; ...
            [data_type_str]})
        
        
         % Save figure  - output -
         fg_name = [data_version ...
             '_' 'PROJ#LMB' '_' datestr(itime,'yyyymmdd') '_SSSerror.' fg_format];
         
         % store figs in this folder
         folder_out_this = [folder_out '/SSSerror'];
         foldercheck_raf(folder_out_this); %! make folder_out
         
         fg_name = [folder_out_this '/' fg_name];
         
         fg_exist = exist(fg_name,'file'); % check fn existence
         
         if save_fig == 1 && fg_exist == 0
             
             save_raf(gcf,fg_name,fg_format); close
             
         end
         
    end
    
    
end


