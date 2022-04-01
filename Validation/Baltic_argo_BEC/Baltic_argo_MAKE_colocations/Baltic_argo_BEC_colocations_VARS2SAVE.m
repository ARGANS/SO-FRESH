% Baltic_argo_BEC_colocations_VARS2SAVE (snippet)
%
% script to pre-locate variables to be colocated by
% Baltic_argo_BEC_colocations.m
%
% version: v1r1 - reduce number of inputs (i.e. colocation one-satellite at
% the time)
%
% History
% v1r0 - Creation of this script
%
% =========================================================================

TT = nan(1,nTOT);
TT2 = nan(ndepth,nTOT);
TT3 = nan(nele,nTOT);

sss_irange = TT3;

sss_error_irange = TT3;

lon_sss_irange = TT3;
lat_sss_irange = TT3;

sss_irange_mn = TT;

sss_irange_md = TT;

sss_error_irange_mn = TT;

sss_error_irange_md = TT;


% Colocation stats (dSSS = Sat-Ref) (X1)
% dSSS_irange_nobs = TT;
% dSSS_irange_mean = TT;
% dSSS_irange_median = TT;
% dSSS_irange_std = TT;
% dSSS_irange_min = TT;
% dSSS_irange_max = TT;
% dSSS_irange_SEM = TT;
% dSSS_irange_Q = nan(nTOT,5);
% dSSS_irange_Q2 = TT;
% dSSS_irange_Q4 = TT;
% dSSS_irange_IQR = TT;
% dSSS_irange_STATS_str = NaN;

% Keep Argo information
% platform identification
platform_irange = TT;
ID_irange = TT;
lon_irange = TT;
lat_irange = TT;
JULD_irange = TT;
count_bin = TT;


% keep ARGO profiles
if bin_data == 0
    % raw profiles
    PRES_irange = TT2;
    SALT_irange = TT2;
    TEMP_irange = TT2;
    
    % interpolated raw profiles
    PRES_intp = TT2;
    SALT_intp = TT2;
    TEMP_intp = TT2;
    
elseif bin_data == 1
    % raw profiles
    PRES_irange = TT;
    SALT_irange = TT;
    TEMP_irange = TT;
    
    % interpolated raw profiles
    PRES_intp = TT;
    SALT_intp = TT;
    TEMP_intp = TT;
    
end