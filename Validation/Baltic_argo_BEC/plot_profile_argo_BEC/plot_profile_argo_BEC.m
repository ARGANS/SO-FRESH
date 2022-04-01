function plot_profile_argo_BEC(TEMP1,SALT1,PRES1,TEMP2,SALT2,PRES2,platform,iTIME,str1,str2,varargin)
%
% Syntax: plot_argo_BEC(TEMP1,SALT1,PRES1,TEMP2,SALT2,PRES2)
%
% Description
% Plot Argo profiles with the RAW and the VERTICAL INTERPOLATED profiles
%
% Input
% TEMP1 = Temperature RAW profile
% SALT1 = Salinity RAW profile
% PRES1 = RAW Presure
% TEMP2 = Temperature vertical interpolated profile
% SALT2 = Salinity vertical interpolated profile
% PRES2 = Presure vertical interpolated
% platform = Argo ID number
% iTIME = time recorded in the file name (see BEC-ARGO files)
% str1 = string for legend first dataset (TEMP1 and SALT1)
% str2 = string for legend second dataset (TEMP2 and SALT2)
% boundedplot_flag (optional):  [1] make boundedplot around mean and std, or do
% not make boundedplot (make individual profiles -- default) [0]
%
% Dimentions of all variables must be: number of depth x number of profiles
%
%
% Output
% N/A
%
% =========================================================================
% Author: rcatany
% current version: v1r0 (2019/11/14)
% =========================================================================






if ~isempty(varargin)
    boundedplot_flag = varargin{1};
else
    boundedplot_flag = 0;
    
end
    


itime_start = min(iTIME(:));
itime_end = max(iTIME(:));

% To check number of profiles (nprof1 and nprof2) in both datasets

if size(TEMP1) == size(SALT1)
    [~,nprof1] = size(TEMP1);
end

if size(TEMP2) == size(SALT2)
    [~,nprof2] = size(TEMP2);
end

if nprof1 == nprof2
    nprof = nprof1;
else
    nprof = [];
end

% ==================
T1 = TEMP1 (:,1:nprof1);
S1 = SALT1 (:,1:nprof1);
P1 = PRES1 (:,1:nprof1);

T2 = TEMP2 (:,1:nprof2);
S2 = SALT2 (:,1:nprof2);
P2 = PRES2 (:,1:nprof2);


plat = unique(platform);

if length(plat)>1
    error('Profiles must be from the same platform number')
end

% multiply the xlim (Tmin,Tmax and Smin,Smax) by percent number (below/above)
percent_number = 0.05; % five percent scalling

factor_lim_min = 1-percent_number; % 0.05 % below the min
factor_lim_max = 1+percent_number; % 0.05 % above the max

Tmin = min([min(T1(:)) min(T2(:))]).*factor_lim_min;
Tmax = max([max(T1(:)) max(T2(:))]).*factor_lim_max;
Smin = min([min(S1(:)) min(S2(:))]).*(factor_lim_min);
Smax = max([max(S1(:)) max(S2(:))]).*(factor_lim_max);
Pmin = min([PRES2(:).*(-1)]);
Pmax = max([PRES2(:).*(-1)]);


% ================
%% Plot ARGO profiles from the same platform
% figure with one profile (case#01) or multiple profiles (case#02)
gcf; hold on;

if boundedplot_flag == 0
    
    marker_sz = 15;
    
    
    % Case#01: figure with one profile
    subplot(1,2,1)
    h1 = plot (T1, -P1,'ko','MarkerFaceColor','r','MarkerSize',marker_sz);
    hold on
    h2 = plot (T2, -P2,'k.','MarkerSize',marker_sz);
    
    xlim([Tmin Tmax])
    ylim([Pmin Pmax])
    
    xlabel ('T (\circ)'); ylabel('depth (m)')
    grid on
    legend ([h1(1);h2(1)],{str1;str2},'location','SouthEast')
    
    if itime_start ~= itime_end
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd') '-' datestr(itime_end,'yyyymmdd')]});
    else
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd')]});
    end
    
    subplot(1,2,2)
    h1 = plot (S1, -P1,'ko','MarkerFaceColor','b','MarkerSize',marker_sz);
    hold on
    h2 = plot (S2, -P2,'k.','MarkerSize',marker_sz);
    
    xlim([Smin Smax])
    ylim([Pmin Pmax])
    
    xlabel ('S (psu)'); ylabel('depth (m)')
    grid on
    legend ([h1(1);h2(1)],{str1; str2},'location','SouthWest')
    if itime_start ~= itime_end
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd') '-' datestr(itime_end,'yyyymmdd')]});
    else
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd')]});
    end
    
    
elseif boundedplot_flag == 1 % if platform multiple profile. Plot average prof.

    % Case#02: BOUNDEDPLOT - figure multiple profiles
    
    T1_mean = nanmean(T1,2);
    S1_mean = nanmean(S1,2);
    
    std_flag = 1; % compute normalized std (n-1); or not (std_flag = 0)
    
    T1_std = nanstd(T1,std_flag,2);
    S1_std = nanstd(S1,std_flag,2);
    
    % value from 0 (trasparent) to 1 (opaque) (type help boundedline)
    shade_transparent = 0.2;
    
    
    subplot(1,2,1)
    h2 = boundedline(T1_mean,-P1,T1_std,...
        '-r','alpha','orientation','horiz',...
        'transparency', shade_transparent); hold on
    
    h1 = plot (T1,-P1,'-','linewidth',2,'Color',[1 1 1].*0.5);
    
    xlim([Tmin Tmax])
    ylim([Pmin Pmax])
    
    xlabel ('T (\circ)'); ylabel('depth (m)')
    grid on
    legend ([h1(1);h2(1)],{str1; str2},'location','SouthWest')
    if itime_start ~= itime_end
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd') '-' datestr(itime_end,'yyyymmdd')]});
    else
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd')]});
    end
    
    % [2] Salinity profile
    subplot(1,2,2)
    h2 = boundedline(S1_mean,-P1,S1_std,...
        '-b','alpha','orientation','horiz',...
        'transparency', shade_transparent); hold on
    
    h1 = plot (S1, -P1,'-','linewidth',2,'Color',[1 1 1].*0.5);
    
    xlim([Smin Smax])
    ylim([Pmin Pmax])
    
    xlabel ('S (psu)'); ylabel('depth (m)')
    grid on
    legend ([h1(1);h2(1)],{str1; str2},'location','SouthWest')
    
    if itime_start ~= itime_end
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd') '-' datestr(itime_end,'yyyymmdd')]});
    else
        title({['platform: ' num2str(plat)];...
            [datestr(itime_start,'yyyymmdd')]});
    end
    
else
    error(['boundedplot_flag must be 0 (plot individual profiles) '...
        'or 1 (plot boundedline, standard deviation around mean)'])
    
    
    
end










