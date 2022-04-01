function [h1,fg_name] = plot_map_model_BEC(lon,lat,param_in,param_str,ibasin,itime)

% Syntax: plot_map_model_BEC(lon,lat,param_in,ibasin)
% (function-snipet of Baltic_model_BEC_analyses.m)
%
% Description
% Plot location of NEMO (model) data within a given region. Model data is
% used as provided by BEC.
%
% Input
% lon/lat
% param_in:
% [1] SSS: Sea Surface Salinity
% [2] SST: Sea Surface Temperature
% [3] MLD: Mixed Layer Depth
% [4] ibasin: Basin number (09: Baltic)
% [5] itime: matlab time number
%
% Output
% h1 is the handle of the figure
% fg_name: figure output filename (string)
%
%
%
% ========================================================================


% Check dimension input variable is equal
[a1,b1,c1] = size(lon);
[a2,b2,c2] = size(lat);
[a3,b3,c3] = size(param_in);

A = ~isequal(a1,a2,a3);
B = ~isequal(b1,b2,b3);
C = ~isequal(c1,c2,c3);

if any (A) || any(B)
    error('lon lat dimensions must be equal')
    
elseif any(C)
    error('param_in must have same dimensions as lon_in and lat_in')
    
end

[xmin,xmax,ymin,ymax,basin_str] = map_lim_raf (ibasin);

Smin = 0;
Smax = 10;

Tmin = -3;
Tmax = 15;

MLDmin = 0;
MLDmax = 50;


% =====================
% Make surface plot

h1 = gcf;
set(gcf,'DefaultAxesFontSize',24);

map_projection = 'lambert';

lat_min = ymin;
lat_max = ymax;
lon_min = xmin;
lon_max = xmax;

lon_step = 10;
lat_step = 5;

fillmap_super(map_projection,...
    lon_min,lon_max,lat_min,...
    lat_max,lon_step,lat_step);

h1 = pcolorm(lat,lon,param_in);  hold on


%% Set colorscale depending on the plotted variable
if strcmpi(param_str,'sss_model')
    colormap('jet');
    caxis ([Smin Smax])
    cb = colorbar;
    ylabel(cb,[param_str])
    
    fg_name = [basin_str '_SSS_NEMO_' datestr(itime,'yyyymmdd')];
    
elseif strcmpi(param_str,'sst_model')
    colormap('jet');
    caxis ([Tmin Tmax])
    cb = colorbar;
    ylabel(cb,[param_str ' (C\circ)'])
    
    fg_name = [basin_str '_SST_NEMO_' datestr(itime,'yyyymmdd')];
    
elseif strcmpi(param_str,'MLD_model')
    colormap('jet');
    caxis ([MLDmin MLDmax])
    cb = colorbar;
    ylabel(cb,[param_str ' (m)'])
    
    fg_name = [basin_str '_MLD_NEMO_' datestr(itime,'yyyymmdd')];
end

