%% ======================================================| DGP1 |==============================================
%% Figure 1: RW+jumps (DGP1)
clear;clc;close all
load("Results DGP1\Monte_Carlo_DGP1_T_50_rho_95.mat");
iter = 400;
fig = figure('Units','inches','Position',[0 0 5.0 3.0],'Color','white');

% Define colors (colorblind-friendly palette)
colors = [0.00, 0.00, 0.00;      % Black for True
          0.00, 0.00, 1.00;      % Strong Blue for TVP-RW
          1.00, 0.00, 0.00];     % Strong Red for AP-VAR

% Define line styles
lineStyles = {'-', '--', ':'};
% Coefficient labels
coeffLabels = {'\beta_1', '\beta_2', '\beta_3', '\beta_4'};
% Number of observations
T = size(BETA,1);   

% Create subplots
for i = 1:4
    ax = subplot(2, 2, i);
    
    % Extract data
    data = [squeeze(BETA(:,i,iter,1)), squeeze(BETA(:,i,iter,2)), squeeze(BETA(:,i,iter,8))];
    
    % Plot with custom colors and styles
    hold on;
    h = zeros(3,1);
    for j = 1:3
        h(j) = plot(1:T, data(:,j), 'LineWidth', 2.0,'Color', colors(j,:), ...
                   'LineStyle', lineStyles{j});
    end
    hold off;
    
    % Formatting
    title(coeffLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time', 'FontSize', 9);
    ylabel('Coefficient', 'FontSize', 9);
    grid on; box on;
    set(ax, 'FontSize', 8, 'LineWidth', 0.8);
    xlim([1 T]);
    
    % Add legend only to first subplot
    if i == 1
        legend(h, {'True', 'TVP', 'AVP'}, ...
               'Location', 'northwest', ...
               'FontSize', 7, ...
               'Box', 'on');
    end
end


set(gcf, 'Color', 'w');
% Save as PDF 
exportgraphics(fig,'Figure_1.pdf','ContentType','vector');

%% ======================================================| DGP2 |==============================================
%% Figure 2: Non-linear TVP regimes thresholds (DGP2)
clear;clc;close all
load("Results DGP2\Monte_Carlo_DGP2_T_100_rho_95.mat")
iter = 650;
fig = figure('Units','inches','Position',[0 0 5.0 3.0],'Color','white');
colors = [0.00, 0.00, 0.00;
          0.00, 0.00, 1.00;
          1.00, 0.00, 0.00];
lineStyles = {'-', '--', ':'};
coeffLabels = {'\beta_1', '\beta_2', '\beta_3', '\beta_4'};

for i = 1:4
    subplot(2, 2, i);
    data = [squeeze(BETA(:,i,iter,1)), ...
            squeeze(BETA(:,i,iter,2)), ...
            squeeze(BETA(:,i,iter,8))];
    
    hold on;
    h = zeros(3,1);
    for j = 1:3
        h(j) = plot(data(:,j), 'LineWidth', 2.0, ...
                   'Color', colors(j,:), ...
                   'LineStyle', lineStyles{j});
    end
    hold off;
    
    title(coeffLabels{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time', 'FontSize', 9);
    ylabel('Coefficient', 'FontSize', 9);
    grid on; box on;
    set(gca, 'FontSize', 8, 'LineWidth', 0.8);
    
    if i == 1
        legend(h, {'True', 'TVP', 'AVP'}, ...
               'Location', 'best', ...
               'FontSize', 7, ...
               'Box', 'on'); 
    end
end
set(gcf, 'Color', 'w');
% Save as  PDF 
exportgraphics(fig,'Figure_2.pdf','ContentType','vector');


