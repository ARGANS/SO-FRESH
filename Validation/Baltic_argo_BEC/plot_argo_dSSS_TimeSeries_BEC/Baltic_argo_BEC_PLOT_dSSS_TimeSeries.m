% Syntax: Baltic_argo_BEC_PLOT_dSSS_TimeSeries.m (script)
%
% Description
% Make Time Series of the difference between Satellite minus Argo satellite
% data (dSSS, known as 'mean bias').
%
%
% Input
% Colocations Argo-2-Satellite as computed by Baltic_argo_BEC_analyses.m
% (script)
%
%
%
% Data source:
% BEC products for validation
%
%
% current version: v1r1 (2020/01/29)
%
% History
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================





%% Regional study (variability within each study region) [as seen in DUM]
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






%% [5.1] SSS Time series


X1 = sort(time_range)';
Y1 = dSSS_gl_irange_mean;
Z1 = dSSS_gl_irange_SEM;

X2 = sort(time_range)';
Y2 = dSSS_nm_irange_mean;
Z2 = dSSS_nm_irange_SEM;

X3 = sort(time_range)';
Y3 = dSSS_ns_irange_mean;
Z3 = dSSS_ns_irange_SEM;


ind1 = isnan(Y1) | Y1<0 | Z1<0 |isnan(Z1);
ind2 = isnan(Y2) | Y2<0 | Z2<0 |isnan(Z2);
ind3 = isnan(Y3) | Y3<0 | Z3<0 |isnan(Z3);

% X1(ind1) = [];
% Y1(ind1) = [];
% Z1(ind1) = [];
% 
% X2(ind2) = [];
% Y2(ind2) = [];
% Z2(ind2) = [];
% 
% X3(ind3) = [];
% Y3(ind3) = [];
% Z3(ind3) = [];


time1 = min(X1(:));
time2 = max(X1(:));

time_str1 = datestr(date1,'yyyymmdd');
time_str2 = datestr(date2,'yyyymmdd');

% ====================================
figure(1); clf
%set(gcf,'DefaultAxesFontSize',24);


% h1 = boundedline(1:length(X1),Y1,Z1,'k-',...
%     'alpha','orientation','vert','transparency',0.3); hold on
% set(h1,'linewidth',2);

h2 = boundedline(X2,Y2,Z2,'b-',...
    'alpha','orientation','vert','transparency',0.3); hold on
set(h2,'linewidth',2);

h3 = boundedline(X3,Y3,Z3,'r-',...
    'alpha','orientation','vert','transparency',0.3); hold on
set(h3,'linewidth',2);

plot(xlim,[0 0],'k--','linewidth',2)

% lg = legend([h1(1);h2(1);h3(1)],...
%     {' SSS-BEC_{global}';' SSS-BEC_{NM}' ; 'SSS_{NS}'},'fontsize',14);

lg = legend([h2(1);h3(1)],...
    {' SSS-BEC_{NM}' ; 'SSS_{NS}'},'fontsize',14);

set(lg,'Location','northwest')

box on
set(gca,'XMinorTick','on','YMinorTick','on')

grid minor
axis tight


set(gca,'xtick',min(X1(:)):15:max(X1(:)))
datetick('x','mm/yy','keepticks')
title({['Average Satellite minus Argo (mean dSSS)'];...
    [basin_str ' ' time_str1 '-' time_str2]},'fontsize',24);
xlabel('time (months)','fontsize',24);
ylabel([ 'dSSS (psu)' ], 'fontsize',24);

fg_name = [basin_str '_MEAN_dSSS_ARGOvsSAT_'...
    time_str1 '_' time_str2];

% Save figure  - output -
folder_this = [folder_figs 'ARGOvsSAT/TIMESERIES/MEAN/dSSS'];
%%
if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name '.' fg_format];

% check fn existence
fg_exist = exist(fg_name,'file');
if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end

