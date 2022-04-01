% Syntax: baltic_plot_argo_map_BEC (script-snipt of Baltic_argo_BEC_analyses.m)
%
% Description
% Plot location of BEC-ferrybox ship based obs within a given region and zoom-in to the
% flaot location.
%
% Version: v1r2
% ========================================================================

ref_str = 'argo'; % input the name of the ref dataset (i.e. argo, ferrybox, or helcom)

%% [1] Baltic+ Plot ALL ref (in situ) location within region
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
 
folder_this = [folder_figs 'MAPS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end


fillmap_baltic


% set title 
if unique(X1.platform_irange) == 1
    
    platform_this = num2str(X1.platform_irange);
    platform_str = num2str(platform_this);
    
    title({[ref_str ' colocations to Baltic ' version_str];...
        [basin_str ' (' region ', ' ...
        num2str(ymin) '-' ...
        num2str(ymax) '\circ N) ' ];...
        ['ID' platform_str,'profiles: ' num2str(nprof)]
        [itime_start_str '-'...
        itime_end_str]});
    
    % Save figure  - output -
    fg_name = [ref_str '_' platform_str '_' itime_start_str '_' itime_end_str ...
        '_MAP_' region '_'...
        version_str '.' fg_format];
    
elseif unique(X1.platform_irange) > 1
    
    title({[ref_str ' colocations to Baltic ' version_str];...
        [basin_str ' (' region ', ' ...
        num2str(ymin) '-' ...
        num2str(ymax) '\circ N) ' ];...
        ['profiles: ' num2str(nprof)]
        [itime_start_str '-'...
        itime_end_str]});
    
    
    % Save figure  - output -
    fg_name = [ref_str '_' itime_start_str '_' itime_end_str ...
        '_MAP_' region '_'...
        version_str '.' fg_format];

    
end

fg_name = [folder_this fg_name];


hold on

h1 = plotm(a,b,...
    'ko','MarkerFaceColor','r','MarkerSize',5); hold on

    
if exist(fg_name,'file') ~= 2
    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
end


% %% [2] Arctic v2.0 Plot ALL Argo float location within region
% a = X2.lat_irange;
% b = X2.lon_irange;
% c = X2.time_irange;
% 
% version_str = X2.version_str;
% 
% 
% nprof = sum(~isnan(a(:)) & ~isnan(b(:)));
% 
% itime_start = c(1);
% itime_end = c(end);
% itime_start_str = datestr(itime_start,'yyyymmdd');
% itime_end_str = datestr(itime_end,'yyyymmdd');
% 
% % ============================================
% 
% % Save figure  - output -
% fg_name = ['argo_' region '_MAP_'...
%     itime_start_str '_' itime_end_str '_' version_str '.' fg_format];
% 
% 
% folder_this = [folder_figs 'MAPS/'];
% 
% if fg_save == 1
%     foldercheck_raf(folder_this); %! make folder_figs
% end
% 
% fg_name = [folder_this fg_name];
% 
% 
% 
%     
% fillmap_arctic    
%     
%     title({['Argo colocations to Arctic ' version_str];...
%         [basin_str ' (' region ', ' ...
%         num2str(ymin) '-' ...
%         num2str(ymax) '\circ N) ' ];...
%         ['profiles: ' num2str(nprof)]
%         [itime_start_str '-'...
%         itime_end_str]});
%     hold on
%     
%     
%     
%     if strcmpi(region,'NA')
%         h1 = plotm(a,b,...
%             'ko','MarkerFaceColor','r','MarkerSize',10); hold on
%         
%     elseif strcmpi(region,'PA')
%         h1 = plotm(a,b,...
%             'ko','MarkerFaceColor','b','MarkerSize',10); hold on
%         
%     elseif strcmpi(region,'ALL')
%         
%         ind1 = b >=-100 & b <=90; % floats in NA
%         ind2 = b >=-180 & b <= -120 | b > 90;% floats in PA
%         
%         nprof_NA = sum(ind1(:));
%         nprof_PA = sum(ind2(:));
%         
%         h1 = plotm(a(ind1),b(ind1),'ko','MarkerFaceColor','r');
%         h2 = plotm(a(ind2),b(ind2),'ko','MarkerFaceColor','b');
%         
%         lg = legend([h1(1),h2(1)],...
%             {['NA (' num2str(nprof_NA) ')'],...
%             ['PA (' num2str(nprof_PA) ')' ]},...
%             'fontsize',18,'location','SouthEast');
%         
%     end
%     
%     
% %     % Save figure  - output -
% %         folder_this = [folder_figs 'MAPS/'];
% %     
% %         if fg_save == 1
% %             foldercheck_raf(folder_this); %! make folder_figs
% %         end
% %     
% %         fg_name = [folder_this fg_name];
%     
% if exist(fg_name,'file') ~= 2    % check fn existence
%     fg_exist = exist(fg_name,'file');
%     if fg_save == 1 && fg_exist == 0
%         save_raf(gcf,fg_name,fg_format); close
%     end
%     
% end
% 
% 
% %% [2] CCI+SSS (v01.07) Plot ALL Argo float location within region
% a = X3.lat_irange;
% b = X3.lon_irange;
% c = X3.time_irange;
% 
% version_str = X3.version_str;
% 
% 
% nprof = sum(~isnan(a(:)) & ~isnan(b(:)));
% 
% itime_start = c(1);
% itime_end = c(end);
% itime_start_str = datestr(itime_start,'yyyymmdd');
% itime_end_str = datestr(itime_end,'yyyymmdd');
% 
% % ============================================
% 
% % Save figure  - output -
% fg_name = ['argo_' region '_MAP_'...
%     itime_start_str '_' itime_end_str '_' version_str '.' fg_format];
% 
% 
% folder_this = [folder_figs 'MAPS/'];
% 
% if fg_save == 1
%     foldercheck_raf(folder_this); %! make folder_figs
% end
% 
% fg_name = [folder_this fg_name];
% 
% 
% 
%     
% fillmap_arctic    
%     
%     title({['Argo colocations to Arctic ' version_str];...
%         [basin_str ' (' region ', ' ...
%         num2str(ymin) '-' ...
%         num2str(ymax) '\circ N) ' ];...
%         ['profiles: ' num2str(nprof)]
%         [itime_start_str '-'...
%         itime_end_str]});
%     hold on
%     
%     
%     
%     if strcmpi(region,'NA')
%         h1 = plotm(a,b,...
%             'ko','MarkerFaceColor','r','MarkerSize',10); hold on
%         
%     elseif strcmpi(region,'PA')
%         h1 = plotm(a,b,...
%             'ko','MarkerFaceColor','b','MarkerSize',10); hold on
%         
%     elseif strcmpi(region,'ALL')
%         
%         ind1 = b >=-100 & b <=90; % floats in NA
%         ind2 = b >=-180 & b <= -120 | b > 90;% floats in PA
%         
%         nprof_NA = sum(ind1(:));
%         nprof_PA = sum(ind2(:));
%         
%         h1 = plotm(a(ind1),b(ind1),'ko','MarkerFaceColor','r');
%         h2 = plotm(a(ind2),b(ind2),'ko','MarkerFaceColor','b');
%         
%         lg = legend([h1(1),h2(1)],...
%             {['NA (' num2str(nprof_NA) ')'],...
%             ['PA (' num2str(nprof_PA) ')' ]},...
%             'fontsize',18,'location','SouthEast');
%         
%     end
%     
%     
% %     % Save figure  - output -
% %         folder_this = [folder_figs 'MAPS/'];
% %     
% %         if fg_save == 1
% %             foldercheck_raf(folder_this); %! make folder_figs
% %         end
% %     
% %         fg_name = [folder_this fg_name];
%     
% if exist(fg_name,'file') ~= 2    % check fn existence
%     fg_exist = exist(fg_name,'file');
%     if fg_save == 1 && fg_exist == 0
%         save_raf(gcf,fg_name,fg_format); close
%     end
%     
% end
% 
% 
% 
% 
% 
