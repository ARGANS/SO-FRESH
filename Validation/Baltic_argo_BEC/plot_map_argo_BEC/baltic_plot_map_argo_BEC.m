% Syntax: plot_map_argo_BEC (script-snipt of Baltic_argo_BEC_analyses.m)
%
% Description
% Plot location of BEC-ARGO floats within a given region and zoom-in to the
% flaot location.
%
% Version: v1r2
% ========================================================================


%% [1] Arctic+ Plot ALL ref (in situ) location within region
a = X1.lat_irange;
b = X1.lon_irange;
c = X1.time_irange;     

version_str = X1.version_str;

nprof = sum(~isnan(a(:)) & ~isnan(b(:)));

itime_start = c(1);
itime_end = c(end);
itime_start_str = datestr(itime_start,'yyyymmdd');
itime_end_str = datestr(itime_end,'yyyymmdd');

% ============================================

% Save figure  - output -
fg_name = ['argo_' region '_MAP_'...
    itime_start_str '_' itime_end_str '_' version_str '.' fg_format];


folder_this = [folder_figs 'MAPS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];



fillmap_arctic

title({['Argo colocations to Arctic ' version_str];...
    [basin_str ' (' region ', ' ...
    num2str(ymin) '-' ...
    num2str(ymax) '\circ N) ' ];...
    ['profiles: ' num2str(nprof)]
    [itime_start_str '-'...
    itime_end_str]});
hold on

    
    
    if strcmpi(region,'NA')
        h1 = plotm(a,b,...
            'ko','MarkerFaceColor','r','MarkerSize',10); hold on
        
    elseif strcmpi(region,'PA')
        h1 = plotm(a,b,...
            'ko','MarkerFaceColor','b','MarkerSize',10); hold on
        
    elseif strcmpi(region,'ALL')
        
        ind1 = b >=-100 & b <=90; % floats in NA
        ind2 = b >=-180 & b <= -120 | b > 90;% floats in PA
        
        nprof_NA = sum(ind1(:));
        nprof_PA = sum(ind2(:));
        
        h1 = plotm(a(ind1),b(ind1),'ko','MarkerFaceColor','r');
        h2 = plotm(a(ind2),b(ind2),'ko','MarkerFaceColor','b');
        
        lg = legend([h1(1),h2(1)],...
            {['NA (' num2str(nprof_NA) ')'],...
            ['PA (' num2str(nprof_PA) ')' ]},...
            'fontsize',18,'location','SouthEast');
        
    end
    
    
%     % Save figure  - output -
%         folder_this = [folder_figs 'MAPS/'];
%     
%         if fg_save == 1
%             foldercheck_raf(folder_this); %! make folder_figs
%         end
%     
%         fg_name = [folder_this fg_name];
    
if exist(fg_name,'file') ~= 2
    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
end


%% [2] Arctic v2.0 Plot ALL Argo float location within region
a = X2.lat_irange;
b = X2.lon_irange;
c = X2.time_irange;

version_str = X2.version_str;


nprof = sum(~isnan(a(:)) & ~isnan(b(:)));

itime_start = c(1);
itime_end = c(end);
itime_start_str = datestr(itime_start,'yyyymmdd');
itime_end_str = datestr(itime_end,'yyyymmdd');

% ============================================

% Save figure  - output -
fg_name = ['argo_' region '_MAP_'...
    itime_start_str '_' itime_end_str '_' version_str '.' fg_format];


folder_this = [folder_figs 'MAPS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];



    
fillmap_arctic    
    
    title({['Argo colocations to Arctic ' version_str];...
        [basin_str ' (' region ', ' ...
        num2str(ymin) '-' ...
        num2str(ymax) '\circ N) ' ];...
        ['profiles: ' num2str(nprof)]
        [itime_start_str '-'...
        itime_end_str]});
    hold on
    
    
    
    if strcmpi(region,'NA')
        h1 = plotm(a,b,...
            'ko','MarkerFaceColor','r','MarkerSize',10); hold on
        
    elseif strcmpi(region,'PA')
        h1 = plotm(a,b,...
            'ko','MarkerFaceColor','b','MarkerSize',10); hold on
        
    elseif strcmpi(region,'ALL')
        
        ind1 = b >=-100 & b <=90; % floats in NA
        ind2 = b >=-180 & b <= -120 | b > 90;% floats in PA
        
        nprof_NA = sum(ind1(:));
        nprof_PA = sum(ind2(:));
        
        h1 = plotm(a(ind1),b(ind1),'ko','MarkerFaceColor','r');
        h2 = plotm(a(ind2),b(ind2),'ko','MarkerFaceColor','b');
        
        lg = legend([h1(1),h2(1)],...
            {['NA (' num2str(nprof_NA) ')'],...
            ['PA (' num2str(nprof_PA) ')' ]},...
            'fontsize',18,'location','SouthEast');
        
    end
    
    
%     % Save figure  - output -
%         folder_this = [folder_figs 'MAPS/'];
%     
%         if fg_save == 1
%             foldercheck_raf(folder_this); %! make folder_figs
%         end
%     
%         fg_name = [folder_this fg_name];
    
if exist(fg_name,'file') ~= 2    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
end


%% [2] CCI+SSS (v01.07) Plot ALL Argo float location within region
a = X3.lat_irange;
b = X3.lon_irange;
c = X3.time_irange;

version_str = X3.version_str;


nprof = sum(~isnan(a(:)) & ~isnan(b(:)));

itime_start = c(1);
itime_end = c(end);
itime_start_str = datestr(itime_start,'yyyymmdd');
itime_end_str = datestr(itime_end,'yyyymmdd');

% ============================================

% Save figure  - output -
fg_name = ['argo_' region '_MAP_'...
    itime_start_str '_' itime_end_str '_' version_str '.' fg_format];


folder_this = [folder_figs 'MAPS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];



    
fillmap_arctic    
    
    title({['Argo colocations to Arctic ' version_str];...
        [basin_str ' (' region ', ' ...
        num2str(ymin) '-' ...
        num2str(ymax) '\circ N) ' ];...
        ['profiles: ' num2str(nprof)]
        [itime_start_str '-'...
        itime_end_str]});
    hold on
    
    
    
    if strcmpi(region,'NA')
        h1 = plotm(a,b,...
            'ko','MarkerFaceColor','r','MarkerSize',10); hold on
        
    elseif strcmpi(region,'PA')
        h1 = plotm(a,b,...
            'ko','MarkerFaceColor','b','MarkerSize',10); hold on
        
    elseif strcmpi(region,'ALL')
        
        ind1 = b >=-100 & b <=90; % floats in NA
        ind2 = b >=-180 & b <= -120 | b > 90;% floats in PA
        
        nprof_NA = sum(ind1(:));
        nprof_PA = sum(ind2(:));
        
        h1 = plotm(a(ind1),b(ind1),'ko','MarkerFaceColor','r');
        h2 = plotm(a(ind2),b(ind2),'ko','MarkerFaceColor','b');
        
        lg = legend([h1(1),h2(1)],...
            {['NA (' num2str(nprof_NA) ')'],...
            ['PA (' num2str(nprof_PA) ')' ]},...
            'fontsize',18,'location','SouthEast');
        
    end
    
    
%     % Save figure  - output -
%         folder_this = [folder_figs 'MAPS/'];
%     
%         if fg_save == 1
%             foldercheck_raf(folder_this); %! make folder_figs
%         end
%     
%         fg_name = [folder_this fg_name];
    
if exist(fg_name,'file') ~= 2    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
end





