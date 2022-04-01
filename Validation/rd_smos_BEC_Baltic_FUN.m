function [data_output] = rd_smos_BEC_Baltic_FUN(itime,ibasin,data_type)
% 
% Syntax: [data_out] = rd_smos_BEC_Baltic_FUN(itime,ibasin,data_type)
%
% Description
%
% Get SSS-SMOS data in the Baltic sea at the whole basin and within each of
% the Study Regions (ref to Baltic+ DUM).
%
% Input
%
% itime: time as in matlab timenum
% ibasin: Baltic, ibaltic = 9;
%
% data_type for Baltic (type help Baltic_plot_fun):
%           [data_type = 4] Baltic BEC NOM product (v1.0)
%           [data_type = 5] Baltic BEC NS product (v1.0)
%           
%[3] Global BEC product (v001) (default)
%           [4] Baltic+ Nominal (NM)

%
% There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
% [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
% [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
% [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
% [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)
%
% 
%
%
% current version: v1r0 (2020/01/30)
%
% History
% v1r0 | 20200130 - script creation
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================

% Get SMOS-data within a given distance from the center if SREG.
r = 200; % km

% Get SSS-SMOS data
[A] = rd_smos_L4_BEC_v1r2(itime,ibasin,data_type);


sss = A.sss;
sss_error = A.sss_error;
lon = A.lon;
lat = A.lat;

%% Baltic-Regional SSS-data
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


%% Get SSS-SMOS data in each Baltic-region
[a,b,c] = size(sss);

% number of study region
Nsreg = length(lon_sreg);

% Pre-locate variables
sss_draft = reshape(sss,a*b,c);
sss_error_draft = reshape(sss_error,a*b,c);

sss_irange = nan([a*b,c, Nsreg]);
sss_error_irange = sss_irange;

lon_irange = nan(a*b,Nsreg);
lat_irange = nan(a*b,Nsreg);


for nn = 1:length(lon_sreg)
    DIST = Distance(lon_sreg(nn),lat_sreg(nn),lon,lat);
    
    irange = find(DIST <= r);
    Ln = length(irange);
    
    lon_irange(1:Ln,nn) = lon(irange);
    lat_irange(1:Ln,nn) = lat(irange);

    sss_irange(1:Ln,:,nn) = sss_draft(irange,:);
    sss_error_irange(1:Ln,:,nn) = sss_error_draft(irange,:);
end

% Return to original shape SSS-inrange variables
[a,b,c] = size(sss);

sss_irange = reshape(sss_irange,a,b,c,Nsreg);
sss_error_irange = reshape(sss_error_irange,a,b,c,Nsreg);

sss_irange = squeeze(sss_irange);
sss_error_irange = squeeze(sss_error_irange);

lon_irange = reshape(lon_irange,a,b,Nsreg);
lat_irange = reshape(lat_irange,a,b,Nsreg);


% ========================================
% Function output (structure)
tt = struct;

tt.lon = lon;
tt.lat = lat;

tt.sss = sss;
tt.sss_error = sss_error;

% Regional Baltic data and statistics

tt.lat_sreg = lat_sreg;
tt.lon_sreg = lon_sreg;

tt.lat_irange = lat_irange;
tt.lon_irange = lon_irange;

tt.SREG = SREG;
tt.SREG_str = SREG_str;
tt.SREG_str2 = SREG_str2;

tt.sss_irange = sss_irange;
tt.sss_error_irange = sss_error_irange;

data_output = tt; clear tt ind

