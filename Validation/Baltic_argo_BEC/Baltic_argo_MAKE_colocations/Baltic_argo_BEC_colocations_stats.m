function [vars_out] = Baltic_argo_BEC_colocations_stats(fn_in)
%
% Syntax (function):
% [vars_out] = Baltic_argo_BEC_colocations_stats(fn_in)
%
% Description
% Compute Statistics of the Argo-2-Satellite colocations
% whithin a searching distance and within a time window (i.e. ±7-days).
%
% The output from this script is to feed into Arctic_argo_BEC_analyses.m
% This script save colocations [Arctic_argo_BEC_colocations_VARS2SAVE.m]
%
% Argo files were prepared by BEC for validation SSS (Baltic+ and Arctic+)
%
%
% Input
% fn_in: filename with the Colocations data (file created with function
% Arctic_argo_BEC_colocations.m)
%
%
% Output
% vars_out is structure and Matlab files (.mat) with colocated profiles and
% with the statistics of dSSS = satellite MINUS reference (eg. ref as
% Argo).
%
%
% current version: v1r0 (2020/03/04) -
%         [1] converted script into function; Originally this code was
%         embeded into Arctic_argo_BEC_colocations
%         [2] split up colocations and stats in
%             two different scripts. This script does colocation (only)
%             to input in stats.
%
% History
%
% =========================================================================
% Author: rcatany
%
% =========================================================================

% load file with Colocations reference-to-satellite
load([fn_in '.mat']);

path_root = ('/Volumes/Rogue/Data/');

ibasin = 9; % Basin number: 7 (Arctic)
[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);

%% Save Colocation output
% folder_out = [path_root...
%     'SSS/' basin_str...
%     '/BEC/Validation/indata/argo/Colocations/' version_str '/'];
% foldercheck_raf(folder_out);

fn_out = [fn_in '_STATS.mat'];

fn_exist = exist(fn_out,'file');

if fn_exist ~= 2
    
    %% [1] dSSS = SATELLITE - ARGO
    depth_ref = 10; % Reference depth (m), by gral. assumtion 10 m S is homogeneous
    
    argo_vars_out = {...
        'platform_irange','ID_irange',...
        'PRES_irange','SALT_irange','TEMP_irange',...
        'lon_irange','lat_irange','JULD_irange','time_irange'...
        'TEMP_intp','SALT_intp','PRES_intp','SALT_ref'};
    
    sat_vars_out = {...
        'sss_irange',...
        'sss_error_irange',...
        'lon_sss_irange','lat_sss_irange','version_str'};
    
    colocated_vars_out = {...
        'sss_irange_mn',...
        'sss_irange_md',...
        'sss_error_irange_mn'};
    
    
    stats_X1_vars_out = {...
        'dSSS_irange_nobs','dSSS_irange','dSSS_irange_mean',...
        'dSSS_irange_median','dSSS_irange_std','dSSS_irange_stdRobust'...
        'dSSS_irange_min','dSSS_irange_max',...
        'dSSS_irange_SEM','dSSS_irange_Q'...
        'dSSS_irange_Q2','dSSS_irange_Q4',...
        'dSSS_irange_IQR','dSSS_irange_STATS_str'};
    
    
    save_vars_out = [argo_vars_out sat_vars_out colocated_vars_out...
        stats_X1_vars_out];
    
    
    % Pre-locate stats variables
    nprof = length(platform_irange);
    
    TT = nan(1,nprof);
    
    % Colocation stats (dSSS = Sat-Ref) (X1)
    dSSS_irange_nobs = TT;
    dSSS_irange = TT;
    dSSS_irange_mean = TT;
    dSSS_irange_median = TT;
    dSSS_irange_std = TT;
    dSSS_irange_stdRobust = TT;
    dSSS_irange_min = TT;
    dSSS_irange_max = TT;
    dSSS_irange_SEM = TT;
    dSSS_irange_Q = nan(nprof,5);
    dSSS_irange_Q2 = TT;
    dSSS_irange_Q4 = TT;
    dSSS_irange_IQR = TT;
    dSSS_irange_STATS_str = [];
    
    clear TT
    
    SALT_ref = nan(1,nprof);
    
    for ncount = 1:length(PRES_intp(1,:))
        SALT_ref_beta = SALT_intp(find(PRES_intp(:,ncount)<=depth_ref),ncount);
        
        ind = ~isnan(SALT_ref_beta);
        SALT_ref_beta = SALT_ref_beta(ind);
        
        if length(SALT_ref_beta) > 1
            SALT_ref_beta = median(SALT_ref_beta);
            
        end
        
        
        if ~isnan(SALT_ref_beta)
            SALT_ref (ncount) = SALT_ref_beta;
            
        else
            
            S_beta = SALT_intp(:,ncount);
            P_beta = PRES_intp(:,ncount);
            
            ind = isnan(S_beta);
            S_beta(ind) = [];
            P_beta(ind) = [];
            
            S_range = range(S_beta);
            
            % keep profile: if [1] there are at least 2 salinity measures; [2]
            % within the first 30 m depth (we want to compare to surface); [3]
            % range bettween salintiies is lower than 0.1 (i.e. about satellite
            % accuracy)
            if length(S_beta) > 1 && S_range <= 0.1 && max(P_beta) < 30
                SALT_ref(ncount) = S_beta(end);
                
            else
                SALT_ref(ncount) = NaN;
                
            end
        end
        
    end; clear *_beta*
    
    
    
    %% Compute dSSS stats
    %% reminder dSSS = Satellite MINUS Argo
    
    
    if ~isempty (SALT_ref)
        sss_irange_beta = sss_irange;
        
        nobs1 = sum(~isnan(sss_irange_beta(:)));
        
        % Match dimensions of satellite data and the ref data (e.g. Argo)
        [row1,col1] = size(sss_irange_beta);
        [row2,col2] = size(SALT_ref);
        
        a = row1/row2;
        b = col1/col2;
        
        SALT_ref_beta = repmat(SALT_ref,a,b);
        dSSS_irange = sss_irange - SALT_ref_beta;
        
        clear *_beta
        
        
        X1 = dSSS_irange;
        
        % Number of Sat-2-Argo colocations
        dSSS_irange_nobs = nobs1;
        
        
        if any(~isnan(X1(:)))
            [X1_mean,X1_median,...
                X1_std,X1_stdRobust,X1_min,...
                X1_max,...
                X1_SEM,X1_Q,...
                X1_Q2,X1_Q4,...
                X1_IQR,X1_STATS_str] =...
                validation_stats(X1);
            
            % Keep STATS output
            % stats output 1- Nominal Sampling
            dSSS_irange_mean = X1_mean;
            dSSS_irange_median = X1_median;
            dSSS_irange_std = X1_std;
            dSSS_irange_stdRobust = X1_stdRobust;
            dSSS_irange_min = X1_min;
            dSSS_irange_max  = X1_max;
            dSSS_irange_SEM  = X1_SEM;
            dSSS_irange_Q = X1_Q;
            dSSS_irange_Q2 = X1_Q2;
            dSSS_irange_Q4 = X1_Q4;
            dSSS_irange_IQR = X1_IQR;
            dSSS_irange_STATS_str = X1_STATS_str;
            
        end
        
        clear X1_*
        
        disp({'New file with Argo-2-satellite colocations '; fn_out})
        save(fn_out,save_vars_out{:})
        
        
        % Setup function output
        TT = struct;
        
        for nn = 1:length(save_vars_out)
            this_param = save_vars_out{nn};
            
            eval(['TT.' this_param '= ' this_param ';']);
            
        end
    end
    
else
    
    TT = load (fn_out);
    
end

vars_out = TT;
clear TT


