function [t_intp,s_intp,p_intp] = zinterp1_ARGO(TEMP,SALT,PRES,z_grid)
% 
% Syntax: [t_intp,s_intp,p_intp] = zinterp1_ARGO(TEMP,SALT,PRES,z_grid)
%
% Description
% Interpolate Z-direction (narrow depth levels): Vertical interpolation 
% (i.e. along depth) of an Argo profile (T and S).
%
% Input
% TEMP = Temperature profile
% SALT = Salinity profile
% PRES = Presure
% z_grid = interpolant (in m) resolution
%
% Output
% t_int = interpolated temperature profile
% s_int = interpolated salinity profile
% p_int = interpolated presure profile
%
% =========================================================================
% Author: rcatany
% current version: v1r0 (2019/11/14) 
% =========================================================================


% Interpolation specs: Grid dimensions
pres1 = 0;              % min pres bin
pres2 = max(PRES(:));      % max depth (m)
grid_p = z_grid;        % interpolant distance (m) between depth levels

int_method = 'linear';              % interpolation method (see doc interp1)
pres_p = [pres1:grid_p:pres2]';     % interpolated presure
zlevels = length(pres_p);           % number of depth levels in the interpolated vars



[i,j] = size(TEMP);

t2 = nan(zlevels,j);
s2 = nan(zlevels,j);
p2 = ones(zlevels,j).*pres_p;


% Make PRES to be same size as SALT (and TEMP)
[a1,b1] = size(SALT);
[a2,b2] = size(PRES);

if b1 ~= b2
    PRES = repmat(PRES,1,b1);
end




nprof1 = length(TEMP(1,:));
nprof2 = length(SALT(1,:));

if nprof1 == nprof2 
    Nobs = nprof1;
end




% interpolate profiles z-axis (depth)
for j = 1 : Nobs
    t1 = TEMP(:,j);
    s1 = SALT(:,j);
    p1 = PRES(:,j); % presure levels are the same for all the profiles
    
    % Remove nan
    ind = isnan(t1) | isnan(p1) | isnan(s1) | p1 > pres2;
    
    p1(ind) = [];
    t1(ind) = [];
    s1(ind) = [];
    
    [p1, index] = unique(p1);
    
    
    dp = diff(p1);
    
    if min(dp) >= 0.01
        t2_interp = interp1(p1,t1(index),pres_p,int_method);
        t2(:,j) = t2_interp;
        
        s2_interp = interp1(p1,s1(index),pres_p,int_method);
        s2 (:,j) = s2_interp;
    end
    % clear vars within loop
    clear t1 s1 p1 ind
end

t_intp = t2;
s_intp = s2;
p_intp = p2;

clear t2 s2 p2