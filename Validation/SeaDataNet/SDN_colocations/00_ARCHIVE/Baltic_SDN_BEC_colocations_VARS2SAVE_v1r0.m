% Baltic_argo_BEC_colocations_VARS2SAVE (snippet)
%
% script to pre-locate variables to be colocated by
% Baltic_argo_BEC_colocations.m
%
% version: v1r0
%
% =========================================================================

TT = nan(1,nTOT);
TT2 = nan(ndepth,nTOT);
TT3 = nan(nele,nTOT);

sss_nm_irange = TT3;
sss_ns_irange = TT3;
sss_gl_irange = TT3;
sss_model_irange = TT3;

sss_nm_error_irange = TT3;
sss_ns_error_irange = TT3;

lon_sss_irange = TT3;
lat_sss_irange = TT3;


sss_nm_irange_mn = TT;
sss_ns_irange_mn = TT;
sss_gl_irange_mn = TT;
sss_model_irange_mn = TT;


sss_nm_irange_md = TT;
sss_ns_irange_md = TT;
sss_gl_irange_md = TT;
sss_model_irange_md = TT;


sss_nm_error_irange_mn = TT;
sss_ns_error_irange_mn = TT;
sss_gl_error_irange_mn = TT;
sss_model_error_irange_mn = TT;


sss_nm_error_irange_md = TT;
sss_ns_error_irange_md = TT;
sss_gl_error_irange_md = TT;
sss_model_error_irange_md = TT;


% Nominal Sampling (X1)
dSSS_nm_irange_nobs = TT;
dSSS_nm_irange_mean = TT;
dSSS_nm_irange_median = TT;
dSSS_nm_irange_std = TT;
dSSS_nm_irange_min = TT;
dSSS_nm_irange_max = TT;
dSSS_nm_irange_SEM = TT;
dSSS_nm_irange_Q = nan(nTOT,5);
dSSS_nm_irange_Q2 = TT;
dSSS_nm_irange_Q4 = TT;
dSSS_nm_irange_IQR = TT;
dSSS_nm_irange_STATS_str = NaN;

% Nodal Sampling (X2)
dSSS_ns_irange_nobs = TT;
dSSS_ns_irange_mean = TT;
dSSS_ns_irange_median = TT;
dSSS_ns_irange_std = TT;
dSSS_ns_irange_min = TT;
dSSS_ns_irange_max = TT;
dSSS_ns_irange_SEM = TT;
dSSS_ns_irange_Q = nan(nTOT,5);
dSSS_ns_irange_Q2 = TT;
dSSS_ns_irange_Q4 = TT;
dSSS_ns_irange_IQR = TT;
dSSS_ns_irange_STATS_str = NaN;


% Global product (X3)
dSSS_gl_irange_nobs = TT;
dSSS_gl_irange_mean = TT;
dSSS_gl_irange_median = TT;
dSSS_gl_irange_std = TT;
dSSS_gl_irange_min = TT;
dSSS_gl_irange_max = TT;
dSSS_gl_irange_SEM = TT;
dSSS_gl_irange_Q = nan(nTOT,5);
dSSS_gl_irange_Q2 = TT;
dSSS_gl_irange_Q4 = TT;
dSSS_gl_irange_IQR = TT;
dSSS_gl_irange_STATS_str = NaN;

% Model product (X4)
dSSS_model_irange_nobs = TT;
dSSS_model_irange_mean = TT;
dSSS_model_irange_median = TT;
dSSS_model_irange_std = TT;
dSSS_model_irange_min = TT;
dSSS_model_irange_max = TT;
dSSS_model_irange_SEM = TT;
dSSS_model_irange_Q = nan(nTOT,5);
dSSS_model_irange_Q2 = TT;
dSSS_model_irange_Q4 = TT;
dSSS_model_irange_IQR = TT;
dSSS_model_irange_STATS_str = NaN;


% Keep SeaDataNet information
lon_irange = TT;
lat_irange = TT;
time_irange = TT;

platform_type_irange = char(TT);
platform_ID_irange = {TT};

% raw profiles
PRES_irange = TT2;
SALT_irange = TT2;
TEMP_irange = TT2;

PRES_intp = TT2;
SALT_intp = TT2;
TEMP_intp = TT2;
