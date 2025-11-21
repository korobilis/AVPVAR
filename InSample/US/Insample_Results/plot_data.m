% Plot Y, Z, and Cumulative Z: U.S. monthly data (Appendix)
clear; clc; close all;
load("insample_results_p1.mat");


figure('Units','inches','Position',[0 0 6.0 3.4],'Color','white'); 
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

% Define crisis periods
financial_crisis_start = datetime('2007-12-01');financial_crisis_end   = datetime('2009-06-30');
covid_start            = datetime('2020-02-01');covid_end              = datetime('2020-06-01');
y_colors = [0 0 0; 0.4 0.4 0.4];   
z_colors = gray(size(Z,2));       

% Line styles
line_styles = {'-','--','-.',':','-'};
A = [1 2 3 5]; % indices of Z variables to plot
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

%% Panel A: VAR Variables
nexttile; hold on;
xregion(financial_crisis_start,financial_crisis_end,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'HandleVisibility','off');
xregion(covid_start,covid_end,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'HandleVisibility','off');

hY(1) = plot(dates,Y(:,1),'Color',y_colors(1,:),'LineWidth',2.2,'LineStyle','-');
hY(2) = plot(dates,Y(:,2),'Color',y_colors(2,:),'LineWidth',2.0,'LineStyle','--');

ylabel('Growth Rate (%)','FontSize',9,'FontWeight','bold');
ylim([-11 9]); grid on;
set(gca,'GridAlpha',0.3,'GridColor',[0.8 0.8 0.8],'FontSize',9,'Box','on','LineWidth',0.8,'XTickLabel',[]);
legend(hY,{'INDPRO','PCEPI'},'Location','southwest','FontSize',8,'Box','on');
title('Panel A: VAR Variables','FontSize',10,'FontWeight','bold');
hold off;

%% Panel B: Uncertainty Indices (Levels)
nexttile; hold on;
xregion(financial_crisis_start,financial_crisis_end,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'HandleVisibility','off');
xregion(covid_start,covid_end,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'HandleVisibility','off');

for i = 1:4
    ii = A(i);
    hZ(i) = plot(dates,Z(:,ii),'Color',z_colors(ii,:),'LineWidth',1.8,'LineStyle',line_styles{ii});
end

ylabel('Z_t','FontSize',9,'FontWeight','bold');
ylim([-7 8]); grid on;
set(gca,'GridAlpha',0.3,'GridColor',[0.8 0.8 0.8],'FontSize',9,'Box','on','LineWidth',0.8,'XTickLabel',[]);
legend(hZ,{'USEPU','GPR','GECON','EBP'},'Location','southwest','FontSize',8,'Box','on');
title('Panel B: Drivers','FontSize',10,'FontWeight','bold');
hold off;

%% Panel C: Cumulative Uncertainty Indices
nexttile; hold on;
xregion(financial_crisis_start,financial_crisis_end,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'HandleVisibility','off');
xregion(covid_start,covid_end,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'HandleVisibility','off');

for i = 1:4
    ii = A(i);
    hZc(i) = plot(dates,cumsum(Z(:,ii)),'Color',z_colors(ii,:),'LineWidth',1.8,'LineStyle',line_styles{ii});
end

xlabel('Time','FontSize',11,'FontWeight','bold');
ylabel('Cumulative Z_t','FontSize',9,'FontWeight','bold');
ylim([-131 90]); grid on;
set(gca,'GridAlpha',0.3,'GridColor',[0.8 0.8 0.8],'FontSize',9,'Box','on','LineWidth',0.8);
xtickformat('yyyy'); xtickangle(0);
legend(hZ,{'USEPU','GPR','GECON','EBP'},'Location','southwest','FontSize',8,'Box','on');
title('Panel C: Cumulative Drivers','FontSize',10,'FontWeight','bold');
hold off;

%% Export for LaTeX
exportgraphics(gcf,'Figures\plots_US.pdf','ContentType','vector','BackgroundColor','none');
