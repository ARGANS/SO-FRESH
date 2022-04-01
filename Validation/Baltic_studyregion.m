function [lon_sreg,lat_sreg,SREG_str,SREG_str2,Marker_Color_sreg] = Baltic_studyregion
%
% Syntax 
% [lat_sreg,lon_sreg,SREG_str,SREG_str2, Marker_Color_sreg] = Baltic_studyregion ()
%
% Description
%
% Load the location (center) of each of the study region in the Baltic
%
% Input 
% N/A
%
% 
% Output
%
% lat_sreg, lon_sreg: center of each study region
% SREG_str: study region string for legend
% SREG_str: study region string for titles
% 
%
% Model Regional study (variability within each study region)
% There are four Baltic+ study regions (see Baltic+ DUM, p. 28):
% [1] Arkona Basin           [ArB] (55°31'30.22"N,  16°16'15.06"E)
% [2] Bothian Sea            [BOS] (61°55'50.19"N,  19°14'46.24"E)
% [3] Gulf of Finland        [GOF] (59°35'55.07"N,  23°07'27.96"E)
% [4] Northern Baltic Proper [NBP] (57°36'10.05"N,  19°48'25.78"E)
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================

% lat and lon as taken from Google Earth


lon_sreg = [...
    16 16 15.06;...
    19 14 46.24;...
    23 07 27.96;...
    19 48 25.78];

lat_sreg = [...
    55 31 30.22;...
    61 55 50.19;...
    59 35 55.07;...
    57 36 10.05];

% convert angles [deg-min-sec] to [deg]
lon_sreg = dms2degrees(lon_sreg);
lat_sreg = dms2degrees(lat_sreg);


% Load Baltic-SMOS grid
%% [1.2.2] Baltic-Nodal Sampling
data_type = 5; % Baltic-Nodal Sampling
itime = datenum(2012,08,01,0,0,0);
ibasin = 9;
[TT] = rd_smos_L4_BEC_v1r3(itime,ibasin,data_type);

lon_sss = TT.lon;
lat_sss = TT.lat;



[lon_grid,lat_grid] = ...
    griddata_raf(lon_sreg,lat_sreg,lon_sss,lat_sss);

lon_sreg = lon_grid; 
lat_sreg = lat_grid;

clear *grid



SREG = [lon_sreg lat_sreg]; % geolocation centre study regions

SREG_str = {'[01] ArB','[02] BOS','[03] GOF','[04] NBP'};
SREG_str2 = {'ArB','BOS','GOF','NBP'}; % use this for figure names

%% Marker_Color 

Marker_Color_sreg = [...
    0.3020    0.7451    0.9333;...
    1         0         0     ;...
    1.0000    1.0000    0.0667;...
    0.7176    0.2745    1.0000];