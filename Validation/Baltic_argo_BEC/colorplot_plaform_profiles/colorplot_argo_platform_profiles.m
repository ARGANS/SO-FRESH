function [h11,h22] = colorplot_argo_platform_profiles(platform,time_ar,T,S,P,varargin)
%
% Syntax plot_platform_profiles(platform,time_ar,T,S,P,T_intp,S_intp,P_intp)
%
% Description
% Plot Argo profiles from the same platform (platform ID#NUMBER)
%
% Input
% platform: is the unique platform number
% time_ar: refers at the time of the cast (Julian date)
% T, S, P: Temperature, Salinity and Presure
%
% Profiles with improved vertical resolution (optional)
% T_intp, S_intp and S_intp: Z-Interpolated profiles 
% 
%
% =========================================================================
%
% rcatany (2019/10/10)
%
% =========================================================================



% set values to the optional variables
if ~isempty(varargin)
    T_intp = varargin{1};
    S_intp = varargin{2};
    P_intp = varargin{3};
else
    T_intp = T;
    S_intp = S;
    P_intp = P;
end


plat = platform;
time_beta = time_ar;

nprof = length(T(1,:)); % number of profiles

time1 = min(time_beta(:));
time2 = max(time_beta(:));

time_str1 = datestr(time1,'yyyymmdd');
time_str2 = datestr(time2,'yyyymmdd');

Tmin = 0;
Tmax = 20;
Smin = nanmin(S(:))*0.95; % lim minium set to 5% below the salinity minimum
Smax = nanmax(S(:))*1.05; % lim maximum set to 5% above the salinity maximum

Pmin = -nanmax(P(:))*1.20;
Pmax = 0;

% =================================
figure; clf
% [1] Temperature profile
subplot(1,2,1)
h11 = plot (T, -P,'k-','linewidth',2);
hold on
h12 = scatter(T_intp(:),-P_intp(:),30,T_intp(:),...
    'fill','MarkerEdgeColor','k');

a = get(gca, 'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)

colormap('jet')

cb = colorbar;

ylabel(cb,'T (\circ)','fontsize',18)
ylim([Pmin Pmax])
xlim ([Tmin Tmax])

xlabel ('T (\circ)','fontsize',24);
ylabel('depth (m)','fontsize',24)
grid on
grid minor

if time1 ~= time2
    title({['platform: ' num2str(plat(1))...
        ' (' num2str(nprof) ' profiles)'];...
        [time_str1 '-' time_str2]},'fontsize',24);
    
else
    title({['platform: ' num2str(plat(1))...
        ' (' num2str(nprof) ' profiles)'];...
        [time_str1]},'fontsize',24);
    
end

% [2] Salinity profile
subplot(1,2,2)
h21 = plot (S, -P,'k-','linewidth',2);
hold on
h22 = scatter(S_intp(:),-P_intp(:),30,S_intp(:),...
    'fill','MarkerEdgeColor','k');

a = get(gca, 'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)


colormap('jet')

cb = colorbar;

ylabel(cb,'S (psu)','fontsize',18)

ylim([Pmin Pmax])
xlim ([Smin Smax])

xlabel ('S (psu)','fontsize',24);
ylabel('depth (m)','fontsize',24)
grid on
grid minor

if time1 ~= time2
title({['platform: ' num2str(plat(1))...
    ' (' num2str(nprof) ' profiles)'];...
    [time_str1 '-' time_str2]},'fontsize',24);

else
title({['platform: ' num2str(plat(1))...
    ' (' num2str(nprof) ' profiles)'];...
    [time_str1]},'fontsize',24);    
    
end



