













% ===============
%% [2.1] boxplots: Satellite data at float location
% shall we use the median of SSS to be compared againts

% string product name
PRODS = {'SSS_{NS}','SSS_{NM}','SSS_{gl}'};


make_plot = 0; % flag to make plot [1], or not [0]
if make_plot == 1
    
    figure(1); clf; hold on;
    set(gcf,'DefaultAxesFontSize',24);
    set(gca,'TickLabelInterpreter','tex');
    
    bplot = boxplot...
        ([sss_ns_irange,sss_nm_irange,sss_gl_irange],...
        'labels',PRODS);
    
    ylim([0 10]);
    
    ylabel('S')
    
    
    h = findobj(gca, 'type', 'text');
    set(h,'FontSize',25);
    set(h,'VerticalAlignment', 'middle');
    set(gca,'FontSize',25);
    
    title({['Satellite SSS around Argo float (plat ID:'...
        num2str(platform(nn))...
        ' , PROF#' num2str(ID(nn)) ')'];...
        ['Distance: ' num2str(r) ' km']});
    
    fg_name = [fn '_SAT-SSS_BOXPLOT_PLAT'...
        num2str(platform(nn))...
        '_PROF#' num2str(ID(nn))...
        '_R' num2str(r) 'KM.'...
        fg_format];
    
    % Save figure  - output -
    folder_this = [folder_figs 'BOXPLOTS/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name];
    
    
    fg_exist = exist(fg_name,'file'); % check fn existence
    
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
    
end

%% *[2.2] Plot: Argo geolocation
% Plots in one map locations of same platform

make_plot = 0; % flag to make plot [1], or not [0]
if make_plot == 1
    
    % snippet to plot map and zoom-in to Argo location
    plot_map_argo_BEC
    
end





%% Example: Interpolated Vs Raw Profile
%% [2.3.1] plot raw Vs interp T and S of a given Profile (nprof)
TEMP2 = t_intp(:,1);
SALT2 = s_intp(:,1);
PRES2 = p_intp(:,1);

str1 = 'RAW';
str2 = 'intp';


make_plot = 0; % flag to make plot [1], or not [0]
if make_plot == 1
    
    figure(4); clf; hold on;
    set(gcf,'DefaultAxesFontSize',18);
    
    plot_profile_argo_BEC(TEMP1,SALT1,PRES1,...
        TEMP2,SALT2,PRES2,platform(nn),...
        iTIME,str1,str2)
    
    % Save figure  - output -
    fg_name = [fn '_RAWvsINTP_PLAT'...
        num2str(platform(nn))...
        '_PROF#' num2str(ID(nn)) ...
        '.'  fg_format];
    
    % Save figure  - output -
    folder_this = [folder_figs 'RAWvsINT/'];
    
    if fg_save == 1
        foldercheck_raf(folder_this); %! make folder_figs
    end
    
    fg_name = [folder_this fg_name];
    
    % check fn existence
    fg_exist = exist(fg_name,'file');
    if fg_save == 1 && fg_exist == 0
        save_raf(gcf,fg_name,fg_format); close
    end
end

%% [2.3.2] More than one profile-platform, make mean and compare all profiles
if nprof > 1
    [t_intp,s_intp,p_intp] =...
        zinterp1_ARGO(TEMP,SALT,PRES,z_grid);
    
    TEMP1 = t_intp;
    SALT1 = s_intp;
    PRES1 = p_intp;
    
    TEMP2 = nanmean(t_intp,2);
    SALT2 = nanmean(s_intp,2);
    PRES2 = p_intp(:,1);
    
    str1 = 'all-profiles';
    str2 = 'mean';
    
    % ===========
    
    make_plot = 0;
    
    if make_plot == 1
        figure(4); clf; hold on
        set(gcf,'DefaultAxesFontSize',18);
        
        plot_profile_argo_BEC(TEMP1,SALT1,PRES1,...
            TEMP2,SALT2,PRES2,platform(nn),...
            iTIME,str1,str2)
        
        % Save figure  - output -
        fg_name = [fn '_BOUNDEDPLOT_PLAT'...
            num2str(platform(nn)) ...
            '.'  fg_format];
        
        % Save figure  - output -
        folder_this = [folder_figs 'BOUNDEDPLOT/'];
        
        if fg_save == 1
            foldercheck_raf(folder_this); %! make folder_figs
        end
        
        fg_name = [folder_this fg_name];
        
        % check fn existence
        fg_exist = exist(fg_name,'file');
        
        if fg_save == 1 && fg_exist == 0
            save_raf(gcf,fg_name,fg_format); close
        end
    end
    
end



%% Plot STATS Argo vs Satellite
STATS_vars = {'mean dSSS_{NS}','median dSSS_{NM}',...
    'mean dSSS_{NM}','median dSSS_{NM}'};

figure

h1 = plot(X1_mean,'ks','MarkerFaceColor','b','MarkerSize',10); hold on
h2 = plot(X1_median,'b-','LineWidth',5);
h3 = plot(X2_mean,'ko','MarkerFaceColor','r','MarkerSize',10); hold on
h4 = plot(X2_median,'r-','LineWidth',5);

plot(xlim, [0 0],'k--','linewidth',2)

lg = legend([h1,h2,h3,h4],STATS_vars{:},'fontsize',18,'location','SouthEast');

ylim([-2 0.5])

grid on
grid minor
box on

set(gca,'FontSize',20)
xlabel ('Colocation numbers','fontsize',20)
ylabel ('dSSS [psu]','fontsize',20)

title({['Salinity Argo (10 m) minus Satellite '],...
    'Baltic ',...
    num2str(iyear)})

% Save figure  - output -
fg_name = [fn '_BOUNDEDPLOT_dSSS' ...
    '.'  fg_format];

% Save figure  - output -
folder_this = [folder_figs 'dSSS/'];

if fg_save == 1
    foldercheck_raf(folder_this); %! make folder_figs
end

fg_name = [folder_this fg_name];

% check fn existence
fg_exist = exist(fg_name,'file');

if fg_save == 1 && fg_exist == 0
    save_raf(gcf,fg_name,fg_format); close
end

