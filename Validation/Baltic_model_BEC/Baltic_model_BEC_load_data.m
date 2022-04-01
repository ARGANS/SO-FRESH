% Baltic_model_BEC_load_data
%
% Description
%
% Snippet to load model data for Baltic_model_BEC_analyses
%
% Inputs
%
% ibasin: basin numbe as Baltic=9
% time_range: study period, depending on model_data range.
%
%  Based on the availability of data, the inputs of the scripts are:
% *[0] time_range = 2011_2016 * ALL data available
% *[1] time_range = 2011_2013 * Use first time range
% *[2] time_range = 2014_2016 * Use second time range
%
%
% =========================================================================
%
% Author: rcatany
%
% =========================================================================


if time_range == 0
    [TT1] = nemo_rd_BEC_Baltic_FUN (2011,0,0,ibasin,grid_output);
    [TT2] = nemo_rd_BEC_Baltic_FUN (2014,0,0,ibasin,grid_output);
    
    SSS_model = cat(3,TT1.SSS_model,TT2.SSS_model);
    SST_model = cat(3,TT1.SST_model,TT2.SST_model);
    MLD_model = cat(3,TT1.MLD_model,TT2.MLD_model);
    
    TIME_model = cat(1,TT1.TIME_model,TT2.TIME_model);
    
    lon = TT1.lon;
    lat = TT1.lat;
    
    lon_irange = TT1.lon_irange;
    lat_irange = TT1.lat_irange;
    
    SSS_model_irange = cat(2,TT1.SSS_irange,TT2.SSS_irange);
    SST_model_irange = cat(2,TT1.SST_irange,TT2.SST_irange);
    MLD_model_irange = cat(2,TT1.MLD_irange,TT2.MLD_irange);    
    
    % Number of observations in the model
    ind_SSS_model_ALL = ~isnan(SSS_model);
    ind_SST_model_ALL = ~isnan(SST_model);
    ind_MLD_model_ALL = ~isnan(MLD_model);
    
    [a,b,c] = size(ind_SSS_model_ALL);
    
    ind_SSS_model_ALL = reshape(ind_SSS_model_ALL,a*b,c);
    ind_SST_model_ALL = reshape(ind_SST_model_ALL,a*b,c);
    ind_MLD_model_ALL = reshape(ind_MLD_model_ALL,a*b,c);
    
    Nobs_SSS_model_ALL = sum(ind_SSS_model_ALL,1);
    Nobs_SST_model_ALL = sum(ind_SST_model_ALL,1);
    Nobs_MLD_model_ALL = sum(ind_MLD_model_ALL,1);
    
    % Number observations in each subregion (irange)
    ind_SSS_model_irange =  ~isnan(SSS_model_irange);
    ind_SST_model_irange =  ~isnan(SST_model_irange);
    ind_MLD_model_irange =  ~isnan(MLD_model_irange);
    
    Nobs_SSS_model_irange = squeeze(sum(ind_SSS_model_irange,1));
    Nobs_SST_model_irange = squeeze(sum(ind_SST_model_irange,1));
    Nobs_MLD_model_irange = squeeze(sum(ind_MLD_model_irange,1));
    
    
    
    % ===========================
    % STATS computation
    SSS_model_mean_sreg = cat(1,TT1.SSS_mean_sreg,TT2.SSS_mean_sreg);
    SST_model_mean_sreg = cat(1,TT1.SST_mean_sreg,TT2.SST_mean_sreg);
    MLD_model_mean_sreg = cat(1,TT1.MLD_mean_sreg,TT2.MLD_mean_sreg);
    
    SSS_model_median_sreg = cat(1,TT1.SSS_median_sreg,TT2.SSS_median_sreg);
    SST_model_median_sreg = cat(1,TT1.SST_median_sreg,TT2.SST_median_sreg);
    MLD_model_median_sreg = cat(1,TT1.MLD_median_sreg,TT2.MLD_median_sreg);
    
    SSS_model_std_sreg = cat(1,TT1.SSS_std_sreg,TT2.SSS_std_sreg);
    SST_model_std_sreg = cat(1,TT1.SST_std_sreg,TT2.SST_std_sreg);
    MLD_model_std_sreg = cat(1,TT1.MLD_std_sreg,TT2.MLD_std_sreg);
    
    SSS_model_SEM_sreg = SSS_model_std_sreg./(Nobs_SSS_model_irange.^0.5);
    SST_model_SEM_sreg = SST_model_std_sreg./(Nobs_SST_model_irange.^0.5);
    MLD_model_SEM_sreg = MLD_model_std_sreg./(Nobs_MLD_model_irange.^0.5);
    
    SSS_model_mean_sreg_ALL = cat(1,TT1.SSS_mean_sreg_ALL,TT2.SSS_mean_sreg_ALL);
    SST_model_mean_sreg_ALL = cat(1,TT1.SST_mean_sreg_ALL,TT2.SST_mean_sreg_ALL);
    MLD_model_mean_sreg_ALL = cat(1,TT1.MLD_mean_sreg_ALL,TT2.MLD_mean_sreg_ALL);
    
    SSS_model_median_sreg_ALL = cat(1,TT1.SSS_median_sreg_ALL,TT2.SSS_median_sreg_ALL);
    SST_model_median_sreg_ALL = cat(1,TT1.SST_median_sreg_ALL,TT2.SST_median_sreg_ALL);
    MLD_model_median_sreg_ALL = cat(1,TT1.MLD_median_sreg_ALL,TT2.MLD_median_sreg_ALL);
    
    SSS_model_std_sreg_ALL = cat(1,TT1.SSS_std_sreg_ALL,TT2.SSS_std_sreg_ALL);
    SST_model_std_sreg_ALL = cat(1,TT1.SST_std_sreg_ALL,TT2.SST_std_sreg_ALL);
    MLD_model_std_sreg_ALL = cat(1,TT1.MLD_std_sreg_ALL,TT2.MLD_std_sreg_ALL);
    
    SSS_model_SEM_sreg_ALL = SSS_model_std_sreg_ALL./ (Nobs_SSS_model_ALL.^0.5);
    SST_model_SEM_sreg_ALL = SST_model_std_sreg_ALL./ (Nobs_SST_model_ALL.^0.5);
    MLD_model_SEM_sreg_ALL = MLD_model_std_sreg_ALL./ (Nobs_MLD_model_ALL.^0.5);
    
    
    
    
    
elseif time_range == 1
    [TT1] = nemo_rd_BEC_Baltic_FUN (2011,0,0,ibasin,grid_output);
    
    SSS_model = TT1.SSS_model;
    SST_model = TT1.SST_model;
    MLD_model = TT1.MLD_model;
    
    TIME_model = TT1.TIME_model;
    
    lon = TT1.lon;
    lat = TT1.lat;
    
    lon_irange = TT1.lon_irange;
    lat_irange = TT1.lat_irange;
    
    SSS_model_irange = TT1.SSS_irange;
    SST_model_irange = TT1.SST_irange;
    MLD_model_irange = TT1.MLD_irange;
    
    % =====================
    % Number of observations in the model
    ind_SSS_model_ALL = ~isnan(SSS_model);
    ind_SST_model_ALL = ~isnan(SST_model);
    ind_MLD_model_ALL = ~isnan(MLD_model);
    
    [a,b,c] = size(ind_SSS_model_ALL);
    
    ind_SSS_model_ALL = reshape(ind_SSS_model_ALL,a*b,c);
    ind_SST_model_ALL = reshape(ind_SST_model_ALL,a*b,c);
    ind_MLD_model_ALL = reshape(ind_MLD_model_ALL,a*b,c);
    
    Nobs_SSS_model_ALL = sum(ind_SSS_model_ALL,1);
    Nobs_SST_model_ALL = sum(ind_SST_model_ALL,1);
    Nobs_MLD_model_ALL = sum(ind_MLD_model_ALL,1);
    
    % Number observations in each subregion (irange)
    ind_SSS_model_irange =  ~isnan(SSS_model_irange);
    ind_SST_model_irange =  ~isnan(SST_model_irange);
    ind_MLD_model_irange =  ~isnan(MLD_model_irange);
    
    Nobs_SSS_model_irange = squeeze(sum(ind_SSS_model_irange,1));
    Nobs_SST_model_irange = squeeze(sum(ind_SST_model_irange,1));
    Nobs_MLD_model_irange = squeeze(sum(ind_MLD_model_irange,1));
    
    
    
    % ===========================
    % STATS computation
    SSS_model_mean_sreg = TT1.SSS_mean_sreg;
    SST_model_mean_sreg = TT1.SST_mean_sreg;
    MLD_model_mean_sreg = TT1.MLD_mean_sreg;
    
    SSS_model_median_sreg = TT1.SSS_median_sreg;
    SST_model_median_sreg = TT1.SST_median_sreg;
    MLD_model_median_sreg = TT1.MLD_median_sreg;
    
    SSS_model_std_sreg = TT1.SSS_std_sreg;
    SST_model_std_sreg = TT1.SST_std_sreg;
    MLD_model_std_sreg = TT1.MLD_std_sreg;
    
    SSS_model_SEM_sreg = SSS_model_std_sreg./(Nobs_SSS_model_irange.^0.5);
    SST_model_SEM_sreg = SST_model_std_sreg./(Nobs_SST_model_irange.^0.5);
    MLD_model_SEM_sreg = MLD_model_std_sreg./(Nobs_MLD_model_irange.^0.5);
    
    
    SSS_model_mean_sreg_ALL = TT1.SSS_mean_sreg_ALL;
    SST_model_mean_sreg_ALL = TT1.SST_mean_sreg_ALL;
    MLD_model_mean_sreg_ALL = TT1.MLD_mean_sreg_ALL;
    
    SSS_model_median_sreg_ALL = TT1.SSS_median_sreg_ALL;
    SST_model_median_sreg_ALL = TT1.SST_median_sreg_ALL;
    MLD_model_median_sreg_ALL = TT1.MLD_median_sreg_ALL;
    
    SSS_model_std_sreg_ALL = TT1.SSS_std_sreg_ALL;
    SST_model_std_sreg_ALL = TT1.SST_std_sreg_ALL;
    MLD_model_std_sreg_ALL = TT1.MLD_std_sreg_ALL;
    
    SSS_model_SEM_sreg_ALL = SSS_model_std_sreg_ALL./ (Nobs_SSS_model_ALL.^0.5);
    SST_model_SEM_sreg_ALL = SST_model_std_sreg_ALL./ (Nobs_SST_model_ALL.^0.5);
    MLD_model_SEM_sreg_ALL = MLD_model_std_sreg_ALL./ (Nobs_MLD_model_ALL.^0.5);
    
    
    
elseif time_range == 2
    [TT1] = nemo_rd_BEC_Baltic_FUN (2014,0,0,ibasin,grid_output);
    
    SSS_model = TT1.SSS_model;
    SST_model = TT1.SST_model;
    MLD_model = TT1.MLD_model;
    
    TIME_model = TT1.TIME_model;
    
    lon = TT1.lon;
    lat = TT1.lat;
    
    lon_irange = TT2.lon_irange;
    lat_irange = TT2.lat_irange;
    
    
    SSS_model_irange = TT2.SSS_irange;
    SST_model_irange = TT2.SST_irange;
    MLD_model_irange = TT2.MLD_irange;
    
    
    % =====================
    % Number of observations in the model
    ind_SSS_model_ALL = ~isnan(SSS_model);
    ind_SST_model_ALL = ~isnan(SST_model);
    ind_MLD_model_ALL = ~isnan(MLD_model);
    
    [a,b,c] = size(ind_SSS_model_ALL);
    
    ind_SSS_model_ALL = reshape(ind_SSS_model_ALL,a*b,c);
    ind_SST_model_ALL = reshape(ind_SST_model_ALL,a*b,c);
    ind_MLD_model_ALL = reshape(ind_MLD_model_ALL,a*b,c);
    
    Nobs_SSS_model_ALL = sum(ind_SSS_model_ALL,1);
    Nobs_SST_model_ALL = sum(ind_SST_model_ALL,1);
    Nobs_MLD_model_ALL = sum(ind_MLD_model_ALL,1);
    
    % Number observations in each subregion (irange)
    ind_SSS_model_irange =  ~isnan(SSS_model_irange);
    ind_SST_model_irange =  ~isnan(SST_model_irange);
    ind_MLD_model_irange =  ~isnan(MLD_model_irange);
    
    Nobs_SSS_model_irange = squeeze(sum(ind_SSS_model_irange,1));
    Nobs_SST_model_irange = squeeze(sum(ind_SST_model_irange,1));
    Nobs_MLD_model_irange = squeeze(sum(ind_MLD_model_irange,1));
    
    % ===========================
    % STATS computation
    
    SSS_model_mean_sreg = TT2.SSS_mean_sreg;
    SST_model_mean_sreg = TT2.SST_mean_sreg;
    MLD_model_mean_sreg = TT2.MLD_mean_sreg;
    
    SSS_model_median_sreg = TT2.SSS_median_sreg;
    SST_model_median_sreg = TT2.SST_median_sreg;
    MLD_model_median_sreg = TT2.MLD_median_sreg;
    
    SSS_model_std_sreg = TT2.SSS_std_sreg;
    SST_model_std_sreg = TT2.SST_std_sreg;
    MLD_model_std_sreg = TT2.MLD_std_sreg;
    
    SSS_model_SEM_sreg = SSS_model_std_sreg./(Nobs_SSS_model_irange.^0.5);
    SST_model_SEM_sreg = SST_model_std_sreg./(Nobs_SST_model_irange.^0.5);
    MLD_model_SEM_sreg = MLD_model_std_sreg./(Nobs_MLD_model_irange.^0.5);
    
    
    
    SSS_model_mean_sreg_ALL = TT2.SSS_mean_sreg_ALL;
    SST_model_mean_sreg_ALL = TT2.SST_mean_sreg_ALL;
    MLD_model_mean_sreg_ALL = TT2.MLD_mean_sreg_ALL;
    
    SSS_model_median_sreg_ALL = TT2.SSS_median_sreg_ALL;
    SST_model_median_sreg_ALL = TT2.SST_median_sreg_ALL;
    MLD_model_median_sreg_ALL = TT2.MLD_median_sreg_ALL;
    
    SSS_model_std_sreg_ALL = TT2.SSS_std_sreg_ALL;
    SST_model_std_sreg_ALL = TT2.SST_std_sreg_ALL;
    MLD_model_std_sreg_ALL = TT2.MLD_std_sreg_ALL;
    
    SSS_model_SEM_sreg_ALL = SSS_model_std_sreg_ALL./ (Nobs_SSS_model_ALL.^0.5);
    SST_model_SEM_sreg_ALL = SST_model_std_sreg_ALL./ (Nobs_SST_model_ALL.^0.5);
    MLD_model_SEM_sreg_ALL = MLD_model_std_sreg_ALL./ (Nobs_MLD_model_ALL.^0.5);
    
    
else
    error('time_range: [0] 2011-2016; [1] 2011-2013; or [2] 2014-2016')
end



