
clear; clc


iyear = 2011;%2011:2013;
imonth = 1%1:12;
REF_str = 'ferrybox';
fg_format = 'png';

ibasin = 9;
bin_sz = 25;
flag_on = 1; % [flag_on= 1] apply flags as (qc = 49 "good value data"); or not [flag_on = 0]
ROUTE = 2;



% ship routes
SHIPS = {'BalticQueen','FinnMaid','Romantika',...
    'SiljaSernade','Transpaper','Victoria'};


% ===================================
%% set flag_str to notice filename
if flag_on == 1
    flag_str = 'FLAG#ON';
elseif flag_on == 0
    flag_str = 'FLAG#NO';
end

%
bin_sz_str = num2str(bin_sz);

[xmin,xmax,ymin,ymax, basin_str] = map_lim_raf(ibasin);
%% folder settings
path_root = ('/Volumes/Rogue/Data/');


for xx = 1: length(ROUTE)
    iships = ROUTE(xx);
    
    if iships == 1
        ships_str = SHIPS{1};
    elseif iships == 2
        ships_str = SHIPS{2};
    elseif iships == 3
        ships_str = SHIPS{3};
    elseif iships == 4
        ships_str = SHIPS{4};
    elseif iships == 5
        ships_str = SHIPS{5};
    elseif iships == 6
        ships_str = SHIPS{6};
    end
    
    
    
    folder_in = [path_root...
        'SSS/' basin_str...
        '/BEC/Validation/indata/' REF_str...
        '/' REF_str '_mat/BINNED/' flag_str '/' ships_str '/monthly/'];
    % =========
    % folder to store figures
    folder_figs = ['/Volumes/Rogue/scratch/Validation/'];
    %
    folder_figs = [folder_figs basin_str '/' REF_str '/BINNED_PLOTS/' flag_str '/'];
    foldercheck_raf(folder_figs);
    
    
    %* loop through time (year, month and days) to Plot BINNED data
    
    
    % set up counter
    CTR = 0;
    
    for yy = 1:length(iyear)
        iYEAR = iyear(yy);
        
        disp(['Processing REF-2-Satellite colocation '...
            basin_str ' ' num2str(iYEAR)]);
        
        for mm = 1:length(imonth)
            iMONTH = imonth(mm);
            
            disp(['year: ' num2str(iYEAR) ' month: ' sprintf('%02.0f',iMONTH)])
            
            % number of days in month
            idays = calendar(iYEAR,iMONTH);
            idays (idays == 0) = [];
            idays = sort(idays);
            
            for dd = 1:length(idays)
                CTR = CTR+1; % loop counter
                
                iDAY = idays(dd);
                
                % central time
                itime_cnt = datenum(iYEAR,iMONTH,iDAY);
                
                % =========================================
                % define time limits in the bin (also filename)
                ndays = 9; % number of days in a baltic+ file
                N = floor(ndays/2);
                itime_start = itime_cnt-N;
                itime_end = itime_cnt+N; clear N
                % time range in each binned file
                itime_bin = itime_start:itime_end;
                
                time_start_str = datestr(itime_start,'yyyymmdd');
                time_end_str = datestr(itime_end,'yyyymmdd');
                
                % filename with BINNED data (as produced with script
                % Baltic_argo_BEC_binning_v1r0.mat)
                fn = (['Baltic_' ships_str '_BINNED_25KM_' ...
                    time_start_str '_'...
                    time_end_str]);
                
                fn_in = [folder_in fn '.mat'];
                fn_exist = exist(fn_in,'file');
                
                if fn_exist == 2
                    
                    % load REF BINNED data
                    load (fn_in)
                    
                    ind = count_bin == 0;
                    count_bin (ind) = NaN;
                    SALT_bin_mean (ind) = NaN;
                    SALT_bin_median (ind) = NaN;
                    SALT_bin_std (ind) = NaN;
                    TEMP_bin_mean (ind) = NaN;
                    TEMP_bin_median (ind) = NaN;
                    TEMP_bin_std (ind) = NaN;
                    
                    nobs_total = nansum(count_bin(:));
                    nobs_total_str = num2str(nobs_total);
                    
                    
                    clear ind
                    
                    ind = ~isnan(SALT_bin_mean) & ~isnan(SALT_bin_mean);
                    
                    num_bins = sum(ind(:));
                    num_bins_str = num2str(num_bins);
                    
                    % keep number of observations and bins (output)
                    NOBS_TOTAL(CTR) = nobs_total;
                    NUM_BINS(CTR) = num_bins;
                    TIME_CNT(CTR) = time_cnt(1);
                    
                    % ========================================================
                    % [1] PLOT COUNTS in BINS
                    folder_this = [folder_figs 'COUNT/'];
                    foldercheck_raf(folder_this);
                    
                    
                    fg_name = [folder_this fn '_COUNT.' fg_format];
                    
                    fg_exist = exist(fg_name,'file');
                    
                    if fg_exist ~= 2
                        figure
                        fillmap_baltic;
                        pcolorm(lat_bin,lon_bin,count_bin); fillmap_baltic;
                        
                        %* plot original REF location
                        plotm(lat_raw,lon_raw,'ro','markersize',5)
                        
                        %* draw grid
                        plotm(lat_bin,lon_bin,'-','color',[1 1 1].*0.5);
                        plotm(lat_bin',lon_bin','-','color',[1 1 1].*0.5);
                        
                        textm(52,15,{['N_{TOTAL} = ' nobs_total_str],...
                            ['BINS_{TOTAL} = ' num_bins_str]},...
                            'BackgroundColor','w','EdgeColor','k','fontsize',20)
                        
                        caxis([0 10])
                        ch = colorbar;
                        ylabel(ch,'number observations','fontsize',18)
                        
                        title({...
                            ['Number of observations in ' ...
                            bin_sz_str ' x ' bin_sz_str ' '...
                            ' boxes ' basin_str];...
                            [time_start_str '-' time_end_str ...
                            '(time_cnt: ' datestr(time_cnt(1),'yyyymmdd') ')' ]
                            ['<plot_Baltic_argo_BEC_binning_v1r0.mat']},...
                            'interpreter','none');
                        
                        ylabel('latitude [ \circ ]','fontsize',12)
                        xlabel('longitude [ \circ ]','fontsize',12)
                        
                        % Save figure
                        save_raf(gcf,fg_name,fg_format);close;
                        
                    end % fg_exist
                    
                    
                    % ========================================================
                    % [2] PLOT SALT mean in BINS
                    folder_this = [folder_figs 'SALT/'];
                    foldercheck_raf(folder_this);
                    
                    fg_name = [folder_this fn '_SALT_mean.' fg_format];
                    
                    fg_exist = exist(fg_name,'file');
                    
                    if fg_exist ~= 2
                        
                        X = lat_bin;
                        Y = lon_bin;
                        Z = SALT_bin_mean;
                        
                        Smin = 1;
                        Smax = 10;
                        
                        figure
                        fillmap_baltic;
                        pcolorm(X,Y,Z); fillmap_baltic;
                        
                        % plot original REF location
                        plotm(lat_raw,lon_raw,'ko',...
                            'markersize',5)
                        
                        % draw grid
                        plotm(lat_bin,lon_bin,'-','color',[1 1 1].*0.5);
                        plotm(lat_bin',lon_bin','-','color',[1 1 1].*0.5);
                        
                        textm(52,15,{['N_{TOTAL} = ' nobs_total_str],...
                            ['BINS_{TOTAL} = ' num_bins_str]},...
                            'BackgroundColor','w','EdgeColor','k','fontsize',20)
                        
                        colormap('jet')
                        caxis([Smin Smax])
                        ch = colorbar;
                        ylabel(ch,'SSS','fontsize',18)
                        
                        title({...
                            ['SALT mean in ' ...
                            bin_sz_str ' x ' bin_sz_str ' '...
                            ' boxes ' basin_str];...
                            [time_start_str '-' time_end_str ...
                            '(time_cnt: ' datestr(time_cnt(1),'yyyymmdd') ')' ]
                            ['<plot_Baltic_argo_BEC_binning_v1r0.mat']},...
                            'interpreter','none');
                        
                        ylabel('latitude [ \circ ]','fontsize',12)
                        xlabel('longitude [ \circ ]','fontsize',12)
                        
                        % Save figure
                        save_raf(gcf,fg_name,fg_format);close;
                        
                    end % fg_exist
                    
                    % ====================================================
                    % [3] PLOT SALT median in BINS
                    folder_this = [folder_figs 'SALT/'];
                    foldercheck_raf(folder_this);
                    
                    
                    fg_name = [folder_this fn '_SALT_median.' fg_format];
                    
                    fg_exist = exist(fg_name,'file');
                    
                    if fg_exist ~= 2
                        
                        X = lat_bin;
                        Y = lon_bin;
                        Z = SALT_bin_median;
                        
                        Smin = 1;
                        Smax = 10;
                        
                        figure
                        fillmap_baltic;
                        pcolorm(X,Y,Z); fillmap_baltic;
                        
                        % plot original REF location
                        plotm(lat_raw,lon_raw,'ko',...
                            'markersize',5)
                        
                        % draw grid
                        plotm(lat_bin,lon_bin,'-','color',[1 1 1].*0.5);
                        plotm(lat_bin',lon_bin','-','color',[1 1 1].*0.5);
                        
                        textm(52,15,{['N_{TOTAL} = ' nobs_total_str],...
                            ['BINS_{TOTAL} = ' num_bins_str]},...
                            'BackgroundColor','w','EdgeColor','k','fontsize',20)
                        
                        colormap('jet')
                        caxis([Smin Smax])
                        ch = colorbar;
                        ylabel(ch,'SSS','fontsize',18)
                        
                        title({...
                            ['SALT median in ' ...
                            bin_sz_str ' x ' bin_sz_str ' '...
                            ' boxes ' basin_str];...
                            [time_start_str '-' time_end_str ...
                            '(time_cnt: ' datestr(time_cnt(1),'yyyymmdd') ')' ]
                            ['<plot_Baltic_argo_BEC_binning_v1r0.mat']},...
                            'interpreter','none');
                        
                        ylabel('latitude [ \circ ]','fontsize',12)
                        xlabel('longitude [ \circ ]','fontsize',12)
                        
                        % Save figure
                        save_raf(gcf,fg_name,fg_format);close;
                        
                    end % fg_exist
                    
                    % ====================================================
                    % [3] PLOT SALT std in BINS
                    folder_this = [folder_figs 'SALT/'];
                    foldercheck_raf(folder_this);
                    
                    
                    fg_name = [folder_this fn '_SALT_median.' fg_format];
                    
                    fg_exist = exist(fg_name,'file');
                    
                    if fg_exist ~= 2
                        
                        X = lat_bin;
                        Y = lon_bin;
                        Z = SALT_bin_std;
                        
                        Smin = 0;
                        Smax = 1;
                        
                        figure
                        fillmap_baltic;
                        pcolorm(X,Y,Z); fillmap_baltic;
                        
                        % plot original REF location
                        plotm(lat_raw,lon_raw,'ko',...
                            'markersize',5)
                        
                        % draw grid
                        plotm(lat_bin,lon_bin,'-','color',[1 1 1].*0.5);
                        plotm(lat_bin',lon_bin','-','color',[1 1 1].*0.5);
                        
                        textm(52,15,{['N_{TOTAL} = ' nobs_total_str],...
                            ['BINS_{TOTAL} = ' num_bins_str]},...
                            'BackgroundColor','w','EdgeColor','k','fontsize',20)
                        
                        colormap('jet')
                        caxis([Smin Smax])
                        ch = colorbar;
                        ylabel(ch,'std SSS','fontsize',18)
                        
                        title({...
                            ['SALT std in ' ...
                            bin_sz_str ' x ' bin_sz_str ' '...
                            ' boxes ' basin_str];...
                            [time_start_str '-' time_end_str ...
                            '(time_cnt: ' datestr(time_cnt(1),'yyyymmdd') ')' ]
                            ['<plot_Baltic_argo_BEC_binning_v1r0.mat']},...
                            'interpreter','none');
                        
                        ylabel('latitude [ \circ ]','fontsize',12)
                        xlabel('longitude [ \circ ]','fontsize',12)
                        
                        % Save figure
                        save_raf(gcf,fg_name,fg_format);close;
                        
                    end % fg_exist
                    
                    
                    
                    
                end %* if filenanme exist, plot binned data
                
                
            end % loop idays
        end % Loop imonth
    end % LOOPS iYARS
end % Loop iSHIPS

