% Univariate UCSV compared to intercepts from: AVP-VAR; TVP-VAR; TVP-VAR_EB

clear; clc; close all;
load("insample_results_p1.mat");
begin = size(beta_save1,2)-size(beta_save2,2)+2;

% Professional figure setup
set(0, 'DefaultAxesFontName', 'Times New Roman');
coeff = 1;
time_axis = dates(begin+p-1:T_full);

% Equation names and panel labels
eq_names = {'RGDP','HICP','STN'};

%% ==================== FIGURE 1: APV-VAR vs TVP-VAR compared to UCSV ====================
figure('Units','inches','Position',[0 0 6.5 2.5],'Color','white');

for eq = 1:3
    subplot(1,3,eq)

    % Calculate series
    series1 = mean(beta_save1(:,begin-1:end,coeff,eq),1)';
    series2 = mean(beta_save2(:,1:end,coeff,eq),1)';
    series3 = mean(beta_save3(:,begin-1:end,coeff,eq),1)';
    series6 = tauhat(begin-1:end-1,eq);

    series1 = zscore(series1);
    series2 = zscore(series2);
    series3 = zscore(series3);
    series6 = zscore(series6);

    hold on;

    % Crisis periods
    xregion(datetime('2008-03-01'), datetime('2009-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');
    xregion(datetime('2020-01-01'), datetime('2020-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');

    % Plot series1, series3, series6
    h = plot(time_axis,[series1,series3,series6],'LineWidth',2);
    set(h(1),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
    set(h(2),'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    set(h(3),'Color',[0.2 0.2 0.2],'LineStyle',':','LineWidth',2);

    grid on;
    set(gca,'GridAlpha',0.3,'FontSize',11,'Box','on','LineWidth',1);
    xtickangle(45);

    % ↓ Reduce only the X‑axis tick font size
    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 12;

    xlabel('Time','FontSize',10,'FontWeight','bold');
    if eq==1
        ylabel('Coefficient','FontSize',10,'FontWeight','bold');
    end

    % Legend only in the second subplot
    if eq==3
        legend({'AVP-VAR','TVP-VAR','UC-SV'}, ...
            'Location','south','FontSize',9,'Box','on');
    end

    title(sprintf('%s Equation',eq_names{eq}),'FontSize',10,'FontWeight','bold');


    hold off;
end

% Save tightly cropped
exportgraphics(gcf,'Figures\intercepts_EA_FB.pdf','ContentType','vector','BackgroundColor','none');


%% ==================== FIGURE 2: APV-VAR vs TVP-VAR-EB compared to UCSV ====================
figure('Units','inches','Position',[0 0 6.5 2.5],'Color','white');

for eq = 1:3
    subplot(1,3,eq)

    % Calculate series
    series1 = mean(beta_save1(:,begin-1:end,coeff,eq),1)';
    series2 = mean(beta_save2(:,1:end,coeff,eq),1)';
    series3 = mean(beta_save3(:,begin-1:end,coeff,eq),1)';
    series6 = tauhat(begin-1:end-1,eq);

    series1 = zscore(series1);
    series2 = zscore(series2);
    series3 = zscore(series3);
    series6 = zscore(series6);

    hold on;

    % Crisis periods
    xregion(datetime('2008-03-01'), datetime('2009-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');
    xregion(datetime('2020-01-01'), datetime('2020-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');

    % Plot series1, series2, series6
    h = plot(time_axis,[series1,series2,series6],'LineWidth',2);
    set(h(1),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
    set(h(2),'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    set(h(3),'Color',[0.2 0.2 0.2],'LineStyle',':','LineWidth',2);

    grid on;
    set(gca,'GridAlpha',0.3,'FontSize',11,'Box','on','LineWidth',1);
    xtickangle(45);

    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 12;

    xlabel('Time','FontSize',10,'FontWeight','bold');
    if eq==1
        ylabel('Coefficient','FontSize',10,'FontWeight','bold');
    end

    % Legend only in the second subplot
    if eq==3
        legend({'AVP-VAR','TVP-VAR-EB','UC-SV'}, ...
            'Location','south','FontSize',9,'Box','on');
    end

    title(sprintf('%s Equation',eq_names{eq}),'FontSize',10,'FontWeight','bold');

 
    hold off;
end

% Save tightly cropped
exportgraphics(gcf,'Figures\intercepts_EA_EB.pdf','ContentType','vector','BackgroundColor','none');





%% ==================== FIGURE 3: APV-VAR vs TVP-VAR compared to UCSV. UNSTANDARDISED intercepts ====================
figure('Units','inches','Position',[0 0 6.5 2.5],'Color','white');

for eq = 1:3
    subplot(1,3,eq)

    % Calculate series
    series1 = mean(beta_save1(:,begin-1:end,coeff,eq),1)';
    series2 = mean(beta_save2(:,1:end,coeff,eq),1)';
    series3 = mean(beta_save3(:,begin-1:end,coeff,eq),1)';
    series6 = tauhat(begin-1:end-1,eq);

    hold on;

    % Crisis periods
    xregion(datetime('2008-03-01'), datetime('2009-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');
    xregion(datetime('2020-01-01'), datetime('2020-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');

    % Plot series1, series3, series6
    h = plot(time_axis,[series1,series3,series6],'LineWidth',2);
    set(h(1),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
    set(h(2),'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    set(h(3),'Color',[0.2 0.2 0.2],'LineStyle',':','LineWidth',2);

    grid on;
    set(gca,'GridAlpha',0.3,'FontSize',11,'Box','on','LineWidth',1);
    xtickangle(45);

    % ↓ Reduce only the X‑axis tick font size
    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 12;

    xlabel('Time','FontSize',10,'FontWeight','bold');
    if eq==1
        ylabel('Coefficient','FontSize',10,'FontWeight','bold');
    end

    % Legend only in the second subplot
    if eq==3
        legend({'AVP-VAR','TVP-VAR','UC-SV'}, ...
            'Location','south','FontSize',9,'Box','on');
    end

    title(sprintf('%s Equation',eq_names{eq}),'FontSize',10,'FontWeight','bold');


    hold off;
end

% Save tightly cropped
exportgraphics(gcf,'Figures\intercepts_EA_FB_unstand.pdf','ContentType','vector','BackgroundColor','none');


%% ==================== FIGURE 4: APV-VAR vs TVP-VAR-EB compared to UCSV. UNSTANDARDISED intercepts ====================
figure('Units','inches','Position',[0 0 6.5 2.5],'Color','white');

for eq = 1:3
    subplot(1,3,eq)

    % Calculate series
    series1 = mean(beta_save1(:,begin-1:end,coeff,eq),1)';
    series2 = mean(beta_save2(:,1:end,coeff,eq),1)';
    series3 = mean(beta_save3(:,begin-1:end,coeff,eq),1)';
    series6 = tauhat(begin-1:end-1,eq);

    series1 = zscore(series1);
    series2 = zscore(series2);
    series3 = zscore(series3);
    series6 = zscore(series6);

    hold on;

    % Crisis periods
    xregion(datetime('2008-03-01'), datetime('2009-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');
    xregion(datetime('2020-01-01'), datetime('2020-06-30'), ...
        'FaceColor',[0.65 0.65 0.65],'FaceAlpha',0.4,'HandleVisibility','off');

    % Plot series1, series2, series6
    h = plot(time_axis,[series1,series2,series6],'LineWidth',2);
    set(h(1),'Color',[0 0 0],'LineStyle','-','LineWidth',2);
    set(h(2),'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
    set(h(3),'Color',[0.2 0.2 0.2],'LineStyle',':','LineWidth',2);

    grid on;
    set(gca,'GridAlpha',0.3,'FontSize',11,'Box','on','LineWidth',1);
    xtickangle(45);

    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 12;

    xlabel('Time','FontSize',10,'FontWeight','bold');
    if eq==1
        ylabel('Coefficient','FontSize',10,'FontWeight','bold');
    end

    % Legend only in the second subplot
    if eq==3
        legend({'AVP-VAR','TVP-VAR-EB','UC-SV'}, ...
            'Location','south','FontSize',9,'Box','on');
    end

    title(sprintf('%s Equation',eq_names{eq}),'FontSize',10,'FontWeight','bold');

 
    hold off;
end

% Save tightly cropped
exportgraphics(gcf,'Figures\intercepts_EA_EB_unstand.pdf','ContentType','vector','BackgroundColor','none');