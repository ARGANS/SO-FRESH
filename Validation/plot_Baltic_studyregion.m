% Description: plot Baltic Subregions study areas
% 
% Project: Baltic+ Salintiy (ESA)
%
% author: rcatany 
% date:   20200812
%



clear all 
close all
clc

% set up folder_figs storage
folder_figs = '/Volumes/Rogue/scratch/Validation/Baltic/';
fg_format = 'png';


% Plot Baltic sub-regional study regions
REGIONS = {'[01] ARB','[02] BOS','[03] GOF','[04] NBP'};

% load satellite grid
load('/Volumes/Rogue/Data/SSS/Baltic/BEC/Baltic_grid/Baltic_plusv1.0_grid.mat');

% =================
% remove data out side Baltic Sea (lon < 14˚E or lat > 66˚N)
ind = lon_sss<= 14 | lat_sss >= 66;

lon_sss(ind) = NaN;
lat_sss(ind) = NaN;

clear ind

% =================
% make ind for each subregion

lon_alfa = lon_sss; lat_alfa = lat_sss;

ind_reg_01 = find(lon_alfa >= 15 & lon_alfa <= 18 &...
    lat_alfa > 53 & lat_alfa <= 58);

ind_reg_02 = find(lon_alfa >= 18 & lon_alfa <= 23 &...
    lat_alfa > 60 & lat_alfa <= 64);

ind_reg_03 = find(lon_alfa >= 23 & lon_alfa <= 26 &...
    lat_alfa > 58 & lat_alfa <= 62);

ind_reg_04 = find(lon_alfa >= 18 & lon_alfa <= 23 &...
    lat_alfa > 56 & lat_alfa <= 59);

%% make plot
figure
fillmap_baltic; hold on
plotm(lat_sss,lon_sss,'k+','markersize',5); 

h1 = plotm(lat_sss(ind_reg_01),lon_sss(ind_reg_01),'ko','MarkerFaceColor','b','markersize',5);
h2 = plotm(lat_sss(ind_reg_02),lon_sss(ind_reg_02),'ko','MarkerFaceColor','k','markersize',5);
h3 = plotm(lat_sss(ind_reg_03),lon_sss(ind_reg_03),'ko','MarkerFaceColor','r','markersize',5);
h4 = plotm(lat_sss(ind_reg_04),lon_sss(ind_reg_04),'ko','MarkerFaceColor','g','markersize',5);

fillmap_baltic;

lg = legend([h1,h2,h3,h4],REGIONS{:},'fontsize',18);

title(['Baltic sub-regional study areas']);

fg_name = ['Baltic_SUBREGIONS'];

folder_this = [folder_figs 'SUBREGIONS/'];

foldercheck_raf(folder_this)

fg_name = [folder_this fg_name '.' fg_format];

% check fn existence
fg_exist = exist(fg_name,'file');

if fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end
