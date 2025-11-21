%% OOS Forecasting results: All Figures Euro Area data 

clear; close all; clc;
for k = 1:4
suffix = sprintf('p%d',k);
loadFile = sprintf("oos_EA_%s.mat", suffix);
saveFile = sprintf("Appendix_Figures\\Appendix_Forecasting_EA_%s.pdf", suffix);
[R5,R6,R7,R8,R9,R10] = Tables_oos_function(loadFile);


modelNames2 = {'AVP-VAR','CP-VAR','CP-VAR SV','TVP-VAR-EB','TVP-VAR','VAR SVO-t','FAVAR','FAVAR SV'};

% Define horizons to plot (in case you want figures with other forecasting
% horizons)
horizons = [1:8];
x_positions = 1:length(horizons);  % Equal spacing for plotting

figure('Units','inches','Position',[0 0 6.5 4.5],'Color','white');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
% Color scheme
colors = [
    0.0, 0.2, 1.0;   % Strong Blue
    0.8, 0.4, 0.1;   % Orange
    0.2, 0.7, 0.3;   % Green
    0.8, 0.2, 0.2;   % Red
    0.6, 0.2, 0.8;   % Purple
    0.4, 0.4, 0.4;   % Gray
    0.9, 0.6, 0.8;   % Pink
    0.0, 0.6, 0.6    % Teal
];
markers = {'o','s','d','^','v','p','h','*'};
titles_text = {
'Panel (a): Ind Prod - MAE','Panel (b): PCEPI - MAE',...
'Panel (c): Ind Prod - QScore90','Panel (d): PCEPI - QScore90',...
'Panel (e): Ind Prod - QScore10','Panel (f): PCEPI - QScore10'
};
ylabels_text = {'MAE Ratio','MAE Ratio','QScore90 Ratio','QScore90 Ratio','QScore10 Ratio','QScore10 Ratio'};
data_matrices = {R5,R6,R7,R8,R9,R10};
%% Create all 6 panels
for panel = 1:6
    ax = nexttile;
    RR = data_matrices{panel};
    h = size(RR,1);
    hold on;
    % Plot all model series
    plot_handles = gobjects(1,8);
    col_idx = [1,2,3,4,5,6,7,8]; % Skip column 5 (benchmark)
    for i = 1:8
        plot_handles(i) = plot(x_positions, RR(horizons,col_idx(i)), ...
            'LineWidth', 2.0, 'MarkerSize', 6, ...
            'Color', colors(i,:), 'Marker', markers{i}, ...
            'MarkerFaceColor', colors(i,:));
    end
    % Reference line
    yline(1.0,'k--','LineWidth',1.5,'HandleVisibility','off');
    % Formatting
    xlim([0.5, length(horizons)+0.5]);
    xticks(x_positions);
    xticklabels(string(horizons));  % Label with actual horizon values
    grid on;
    set(ax,'GridAlpha',0.25,'GridLineStyle','-','GridColor',[0.8 0.8 0.8]);
    set(ax,'FontSize',10,'Box','on','LineWidth',1.0);
    % Labels
    if panel >= 5
        xlabel('Horizon','FontSize',10,'FontWeight','bold');
    end
    if mod(panel,2)==1
        ylabel(ylabels_text{panel},'FontSize',10,'FontWeight','bold');
    end

    title(titles_text{panel},'FontSize',10,'FontWeight','bold');
    hold off;
end

lgd = legend(modelNames2,'Orientation','horizontal','FontSize',9,'NumColumns',4);
lgd.Layout.Tile = 'south';

%% Save  PDF  
exportgraphics(gcf, saveFile, ...
 'ContentType','vector','BackgroundColor','none');


end