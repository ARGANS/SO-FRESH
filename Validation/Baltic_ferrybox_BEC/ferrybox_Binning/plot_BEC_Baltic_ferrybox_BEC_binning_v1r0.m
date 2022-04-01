
clear; clc


iyear = 2011:2013;
imonth = 1:12;
REF_str = 'ferrybox';
fg_format = 'png';

ibasin = 9;
bin_sz = 25;
flag_on = 1; % [flag_on= 1] apply flags as (qc = 49 "good value data"); or not [flag_on = 0]
ROUTE = [2,5]; % only: 2 or 5
BEC_product = 1; % Nominal [1], or Nodal Sampling [2]


if BEC_product == 1
    product_str = 'NOM';
    
elseif BEC_product == 2
    product_str = 'NS';
    
end


% ship routes
SHIPS = {'BalticQueen','Finnmaid','Romantika',...
    'SiljaSernade','Transpaper','Victoria'};


% ===================================
%% set flag_str to notice filename
if flag_on == 1
    flag_str = 'BEC';
    
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
        '/Colocations/'  flag_str '/'];
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
        
        %* BEC Colocation filename
        fn = (['COLLOCATION_SMOS_'...
            product_str '_' ships_str '_' num2str(iYEAR)]);
        
        fn_in = [folder_in fn '.mat'];
        fn_exist = exist(fn_in,'file');
        
        % load BEC colocations
        
        if exist(fn_in,'file')
            TT = load (fn_in);
            
            itime_BEC = TT.TIME;
            
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
                    
                    ind_time = find(itime_cnt == itime_BEC);
                    
                    if ~isempty(ind_time)
                        lon = TT.LON;
                        lat = TT.LAT;
                        SSS = TT.SSS(:,:,ind_time);
                        TSG = TT.TSG(:,:,ind_time);
                        dSSS = TT.ANOM(:,:,ind_time);
                        
                        SSS (SSS == -999) = NaN;
                        
                        % Make plots (if there is TSG data)
                        
                        ind = ~isnan(TSG);
                        
                        num_bins = sum(ind(:));
                        num_bins_str = num2str(num_bins);
                        
                        % keep number of observations and bins (output)
                        % NOBS_TOTAL(CTR) = nobs_total;
                        NUM_BINS(CTR) = num_bins;
                        TIME_CNT(CTR) = itime_cnt(1);
                        

                        if num_bins > 0
                            
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
                            
                            
                            % ========================================================
                            % [1] PLOT SALT mean in BINS
                            folder_this = [folder_figs 'SALT/'];
                            foldercheck_raf(folder_this);
                            
                            fn_fg = (['COLLOCATION_SMOS_SALT_'...
                                product_str '_' ships_str '_'...
                                time_start_str '_' time_end_str]);
                            
                            fg_name = [folder_this fn_fg '.' fg_format];
                            
                            fg_exist = exist(fg_name,'file');
                            
                            
                           
                            if fg_exist ~= 2
                                
                                X = lat;
                                Y = lon;
                                Z = TSG;
                                
                                Smin = 1;
                                Smax = 10;
                                
                                figure
                                fillmap_baltic;
                                pcolorm(X,Y,Z); fillmap_baltic;
                                
                                %* plot original REF location <BEC
                                %* do not include the raw locations of REF
                                % plotm(lat,lon,'ko','markersize',5)
                                
                                % draw grid
                                plotm(lat,lon,'-','color',[1 1 1].*0.5);
                                plotm(lat',lon','-','color',[1 1 1].*0.5);
                                
                                textm(52,15,{['N_{TOTAL} = < >'],...
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
                                    ' (time_cnt: ' datestr(itime_cnt(1),'yyyymmdd') ')' ]
                                    ['<plot_BEC_Baltic_ferrybox_BEC_binning_v1r0.m']},...
                                    'interpreter','none');
                                
                                ylabel('latitude [ \circ ]','fontsize',12)
                                xlabel('longitude [ \circ ]','fontsize',12)
                                
                                % Save figure
                                save_raf(gcf,fg_name,fg_format);close;
                                
                            end % fg_exist
                            
                            % ========================================================
                            % [2] PLOT SALT mean in BINS
                            folder_this = [folder_figs 'dSSS/'];
                            foldercheck_raf(folder_this);
                            
                            fn_fg = (['COLLOCATION_SMOS_dSSS_'...
                                product_str '_' ships_str '_'...
                                time_start_str '_' time_end_str]);
                            
                            fg_name = [folder_this fn_fg '.' fg_format];
                            
                            fg_exist = exist(fg_name,'file');
                            
                            if fg_exist ~= 2
                                
                                X = lat;
                                Y = lon;
                                Z = dSSS;
                                
                                Smin = -1;
                                Smax = 1;
                                
                                figure
                                fillmap_baltic;
                                pcolorm(X,Y,Z); fillmap_baltic;
                                
                                %* plot original REF location <BEC
                                %* do not include the raw locations of REF
                                % plotm(lat,lon,'ko','markersize',5)
                                
                                % draw grid
                                plotm(lat,lon,'-','color',[1 1 1].*0.5);
                                plotm(lat',lon','-','color',[1 1 1].*0.5);
                                
                                textm(52,15,{['N_{TOTAL} = < >'],...
                                    ['BINS_{TOTAL} = ' num_bins_str]},...
                                    'BackgroundColor','w','EdgeColor','k','fontsize',20)
                                
                                colormap('jet')
                                caxis([Smin Smax])
                                ch = colorbar;
                                colormap(french2);
                                ylabel(ch,'SSS','fontsize',18)
                                
                                title({...
                                    ['SALT mean in ' ...
                                    bin_sz_str ' x ' bin_sz_str ' '...
                                    ' boxes ' basin_str];...
                                    [time_start_str '-' time_end_str ...
                                    ' (time_cnt: ' datestr(itime_cnt(1),'yyyymmdd') ')' ]
                                    ['<plot_BEC_Baltic_ferrybox_BEC_binning_v1r0.m']},...
                                    'interpreter','none');
                                
                                ylabel('latitude [ \circ ]','fontsize',12)
                                xlabel('longitude [ \circ ]','fontsize',12)
                                
                                % Save figure
                                save_raf(gcf,fg_name,fg_format);close;
                                
                            end % fg_exist                            
                            
                            
                            % ========================================================
                            % [3] PLOT SALT mean and SSS in BINS
                            folder_this = [folder_figs 'dSSS2/'];
                            foldercheck_raf(folder_this);
                            
                            fn_fg = (['COLLOCATION_SMOS_dSSS2_'...
                                product_str '_' ships_str '_'...
                                time_start_str '_' time_end_str]);
                            
                            fg_name = [folder_this fn_fg '.' fg_format];
                            
                            fg_exist = exist(fg_name,'file');
                            
                            if fg_exist ~= 2
                                
                                X = lat;
                                Y = lon;
                                z1 = SSS;
                                z2 = TSG;
                                
                                Smin = 1;
                                Smax = 10;
                                
                                figure
                                fillmap_baltic;
                                pcolorm(X,Y,z1); fillmap_baltic;
                                
                                ind2 = ~isnan(z2);
                                scatterm(X(ind2),Y(ind2),50,z2(ind2),...
                                    'filled','MarkerEdgeColor','k')
                                
                                
                                
                                %* plot original REF location <BEC
                                %* do not include the raw locations of REF
                                % plotm(lat,lon,'ko','markersize',5)
                                
                                % draw grid
                                plotm(lat,lon,'-','color',[1 1 1].*0.5);
                                plotm(lat',lon','-','color',[1 1 1].*0.5);
                                
                                textm(52,15,{['N_{TOTAL} = < >'],...
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
                                    ' (time_cnt: ' datestr(itime_cnt(1),'yyyymmdd') ')' ]
                                    ['<plot_BEC_Baltic_ferrybox_BEC_binning_v1r0.m']},...
                                    'interpreter','none');
                                
                                ylabel('latitude [ \circ ]','fontsize',12)
                                xlabel('longitude [ \circ ]','fontsize',12)
                                
                                % Save figure
                                save_raf(gcf,fg_name,fg_format);close;
                                
                            end % fg_exist                            
                                                        
                            
                            
                        end % num_bis>0 (make plot if there is any TSG data)
                        
                        
                    end
                    
                    
                    % ========================================================
                    % [1] PLOT COUNTS in BINS
                    
                    % =================================================
                    % [2] PLOT SALT mean in BINS
                    
                    % =================================================
                    % [3] PLOT SALT median in BINS
                    
                    % =================================================
                    % [4] PLOT SALT std in BINS
                    
                    
                    
                end % loop idays
            end % Loop imonth
            
        end
        
    end % LOOPS iYARS
end % Loop iSHIPS

