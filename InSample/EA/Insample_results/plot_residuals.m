%% Figure In-Sample cumsum(residuals)
clear;clc;close all
load("insample_results_p2.mat");
begin = 42; % TVP-VAR-EB use 40 observations for pre estimation
figure('Units','inches','Position',[0 0 5.5 2.5],'Color','white'); 
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
colors = [0.0, 0.0, 0.0;       % Black for AP-BVAR
          0.4, 0.4, 0.4;       % Dark gray  
          0.7, 0.7, 0.7];      % Light gray 
line_styles = {'-', '--', '-.'};
line_widths = [5.0, 3.5, 3.5];
% Variable names
var_names = {'Real GDP Growth', 'HICP'};
% Crisis periods
if isdatetime(dates)
    gfc_start = datetime('2008-03-01');
    gfc_end = datetime('2009-06-30');
    covid_start = datetime('2020-01-01');
    covid_end = datetime('2020-06-30');
else
    gfc_start = datetime('2008-03-01');
    gfc_end = datenum('2009-06-30');
    covid_start = datenum('2020-01-01');
    covid_end = datenum('2020-06-30');
end
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_counter = 1;

%% Panels: Cumulative Residuals
for kk = 1:2
    nexttile; hold on;
    % Data
    data_matrix = [cumsum(mean(abs(residuals_save1(:,begin-1:end,kk)),1)'),...
                   cumsum(mean(abs(residuals_save3(:,begin-1:end,kk)),1)'), ...      
                   cumsum(mean(abs(residuals_save2(:,1:end,kk)),1)')];
    % Y-limits
    y_max = max(data_matrix(:)) * 1.08; y_min = 0;
    % Crisis shading
    fill([gfc_start gfc_end gfc_end gfc_start],[y_min y_min y_max y_max],...
         [0.85 0.85 0.85],'FaceAlpha',0.4,'EdgeColor','none','HandleVisibility','off');
    fill([covid_start covid_end covid_end covid_start],[y_min y_min y_max y_max],...
         [0.85 0.85 0.85],'FaceAlpha',0.4,'EdgeColor','none','HandleVisibility','off');
    % Plot lines with dates on X-axis
    h = plot(dates(T_full-size(data_matrix,1)+1:T_full), data_matrix, 'LineWidth', 2);
    for i = 1:length(h)
        set(h(i),'Color',colors(i,:),'LineStyle',line_styles{i},'LineWidth',line_widths(i));
    end
    % Legend
    if kk==1
        legend(h,{'AVP-VAR','TVP-VAR','TVP-VAR-EB'},...
            'Location','northwest','FontSize',8,'Box','on','EdgeColor',[0.5 0.5 0.5]);
    end
    % Labels
    if kk==1
        ylabel('Cumulative |Residuals|','FontSize',10,'FontWeight','normal');
    end
    title(var_names{kk},'FontSize',10,'FontWeight','bold');

    grid on;
    set(gca,'GridAlpha',0.25,'GridLineStyle','-','GridColor',[0.8 0.8 0.8]);
    set(gca,'FontSize',9,'LineWidth',0.8,'Box','on');
    xlim([dates(begin+1), dates(T_full)]);
    ylim([y_min, y_max]);
    xtickformat('yyyy'); 
    xtickangle(45);
    hold off;
end
%% Export PDF 
exportgraphics(gcf,'Figures\residuals_EA.pdf','ContentType','vector','BackgroundColor','none');




