
clear;
clc


path_root = '/Volumes/Rogue/';

folder_tmp = ['~/Documents/GitHub/tmp/'];
foldercheck_raf(folder_tmp);

folder_out = folder_tmp; % keep folder out as tmp until new order is given


% inputs script
iyear = 2012;
imonth = 1:12;

ibasin = 9;

[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);


nyear = length(iyear);

% pre-alocate variables of size (a,b):
% a = float number and b = 1;
a = 1;
b = 2000; % number of floats
c= 1200; % number of depth levels

lon_in = nan(a,b);
lat_in = nan(a,b);
plat_in  = nan(a,b);
temp_in = nan(c,b);
salt_in = nan(c,b);
pres_in = nan(c,b);

temp_qc_in = nan(c,b);
salt_qc_in = nan(c,b);
pres_qc_in = nan(c,b);

TEMP_ADJUSTED_ERROR_in = nan(c,b);
SALT_ADJUSTED_ERROR_in = nan(c,b);
PRES_ADJUSTED_ERROR_in = nan(c,b);

vars_save = {'lon_in', 'lat_in','plat_in','temp_in','salt_in','pres_in',...
    'temp_qc_in','salt_qc_in','pres_qc_in','TEMP_ADJUSTED_ERROR_in',...
    'SALT_ADJUSTED_ERROR_in','PRES_ADJUSTED_ERROR_in','xmin','xmax',...
    'ymin','ymax','basin_str','Nobs'};

% file nanme out with the Argo floats in the region
fn_out = [folder_out basin_str '_ARGO_tmp_' num2str(iyear(1))...
    '_' num2str(iyear(end)) '.mat'];

fn_exist = exist(fn_out,'file');

if fn_exist ~= 2
ind_count = 0; % set to zero flout counter within the region
    
    for yy = 1: nyear
        for mm = 1:12
            year_str = num2str(iyear(yy));
            month_str = sprintf('%02.0f',imonth(mm));
            
            % number of days each month
            num_days = calendar(iyear(yy),imonth(mm))';
            
            num_days(num_days == 0) = [];
            num_days = length(num_days);
            
            
            folder_in = [path_root 'Data/Argo/argo_mat/ATL/'...
                year_str '/' month_str '/'];
            
            
            for dd = 1:num_days
                
                fn = [folder_in ...
                    'ATL' year_str month_str sprintf('%02.0f',dd) '.mat'];
                
                load (fn);
                
                
                % Float in range within the Baltic sea (lon < 10 E and lat < 65)
                irange = find (lon >= 10 & lon <= xmax & lat >= ymin & lat <= 65);
                
                if ~isempty (irange)
                    ind_count = ind_count + length(irange); % float counter
                    
                    % number of depth levels and number of floats
                    [ndepth,nfloats] = size(temp(:,irange));
                    
                    if nfloats ~= length(irange)
                        error ('stop')
                    end
                    
                    % index to fill in matrix (ind2)
                    ind = find (isnan(lon_in),1,'first');
                    
                    lon_in (ind : ind-1 + nfloats) = lon(irange);
                    lat_in (ind : ind-1 + nfloats) = lat(irange);
                    plat_in (ind : ind-1 + nfloats) = platform(irange);
                    
                    temp_in (1:ndepth,ind : ind-1 + nfloats) = temp(:,irange);
                    salt_in (1:ndepth,ind : ind-1 + nfloats) = salt(:,irange);
                    pres_in (1:ndepth,ind : ind-1 + nfloats) = pres(:,irange);
                    
                    temp_qc_in (1:ndepth,ind : ind-1 + nfloats) = temp_qc(:,irange);
                    salt_qc_in (1:ndepth,ind : ind-1 + nfloats) = salt_qc(:,irange);
                    pres_qc_in (1:ndepth,ind : ind-1 + nfloats) = pres_qc(:,irange);
                    
                    TEMP_ADJUSTED_ERROR_in (1:ndepth,ind : ind-1 + nfloats) = TEMP_ADJUSTED_ERROR(:,irange);
                    SALT_ADJUSTED_ERROR_in (1:ndepth,ind : ind-1 + nfloats) = PSAL_ADJUSTED_ERROR(:,irange);
                    PRES_ADJUSTED_ERROR_in (1:ndepth,ind : ind-1 + nfloats) = PRES_ADJUSTED_ERROR(:,irange);
                end
            end
        end
    end
    % save file with Argo floats operating within the region
    Nobs = ind_count;   % number of observations
    save(fn_out,vars_save{:});
    
else
    load (fn_out);
end





%% Plot ARGO floats within Baltic region
if Nobs ~= 0
    figure (1); clf; hold on
    
    lat_min2 = ymin;
    lat_max2 = ymax;
    lon_min2 = xmin;
    lon_max2 = xmax;
    
    lat_step = 5;
    lon_step = 10;
    
    axesm('mapprojection','lambertstd','maplatlimit',[lat_min2 lat_max2],'maplonlimit',[lon_min2 lon_max2],...
        'parallellabel','on','plabellocation',[lat_min2:lat_step:lat_max2],'meridianlabel','on',...
        'mlabellocation',[lon_min2:lon_step:lon_max2],'fontsize',14)
    
    framem on
    axis off
    tightmap
    geoshow('landareas.shp','facecolor',[1 1 1]*.5)
    gridm on
    title({[basin_str ' (' num2str(lat_min2) '-' ...
        num2str(lat_max2) '\circ N) ' num2str(iyear(1)) '-' num2str(iyear(end))]; ...
        ['N = ' num2str(Nobs)]})
    hold on
    
    plotm(lat_in,lon_in,'ko','MarkerFaceColor','r')
    
end

%% Interpolate Z-direction (narrow depth levels)
if Nobs ~= 0

% Interpolation specs: Grid dimensions
pres1 = 0;              % min pres bin
pres2 = 500;            % max depth (m)
grid_p = 0.5;             % interpolant distance (m) between depth levels

int_method = 'linear';              % interpolation method (see doc interp1)
pres_p = [pres1:grid_p:pres2]';     % interpolated presure
zlevels = length(pres_p);           % number of depth levels in the interpolated vars

[i,j] = size(temp_in);

t2 = nan(zlevels,j);
s2 = nan(zlevels,j);
p2 = ones(zlevels,j).*pres_p;


% interpolate profiles z-axis (depth)
for j = 1 : Nobs
    t1 = temp_in(:,j);
    s1 = salt_in(:,j);
    p1 = pres_in(:,j); % presure levels are the same for all the profiles
    
    ind = isnan(t1) & isnan(p1) & isnan(s1) | p1 > pres2;
    
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


% ==============================
% Example: results of interpolation
% plot floats in a given location
% locations as:
% #River Mouths
% LOC#01 : Vistula River
% LOC#02 : Neva River
% # Midle Section
% LOC#03 : Gotland to Latvia
% # Exit Section
% LOC#04 
% #North Section
% LOC05 : Gulf of Bothnia
% LOC06 : Gulf of Finland

ind_01 = find(plat_in == 3901940);

lon_01 = lon_in(1,ind_01);
lat_01 = lat_in(1,ind_01);

temp_in01 = temp_in(:,ind_01);
salt_in01 = salt_in(:,ind_01);
pres_in01 = pres_in(:,ind_01);

temp_in_diff = diff(temp_in(:,ind_01));


t_intp01 = t_intp(:,ind_01);
s_intp01 = s_intp(:,ind_01);
p_intp01 = p_intp(:,ind_01);

figure (2); clf; hold on
plot(lon_01,lat_01,'or','markersize',10)
xlim([xmin xmax]); ylim([ymin ymax])
fillmap([0.5 0.5 0.5],1)

figure(3); clf; 
subplot(1,2,1);
h1 = plot (temp_in01,-pres_in01,'ko'); hold on
h2 = plot (t_intp01,-p_intp01,'r.');
xlabel('Temperature (?C)');
ylabel('Depth (m)')

% legend([h1(1),h2(1)],{'raw','interp'},'location','southeast','fontsize',14)

title(['Temperature ARGO profile with z-axis ' ...
    int_method ' interpolation'])

subplot(1,2,2);
h1 = plot (salt_in01,-pres_in01,'ko'); hold on
h2 = plot (s_intp01,-p_intp01,'b.');
xlabel('Salinity (psu)');
ylabel('Depth (m)')

legend([h1(1),h2(1)],{'raw','interp'},'location','southeast','fontsize',14)

title(['Salinity ARGO profile with z-axis ' ...
    int_method ' interpolation'])

end
