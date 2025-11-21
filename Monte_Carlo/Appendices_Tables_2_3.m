%% Monte Carlo results: Tables in Appendices

%% ========================================== |Results Table 2 Appendices (DGP1)| =================================================
clear;clc;
disp('=============================================================== ');
disp('Appendices Table 2 (DGP1) ');
disp('=============================================================== ');
% ------------ |T = 50; rho = 0.00|---------
load("Results DGP1\Monte_Carlo_DGP1_T_50_rho_0.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel A: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 100; rho = 0.00|---------
load("Results DGP1\Monte_Carlo_DGP1_T_100_rho_0.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel B: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 200; rho = 0.00|---------
load("Results DGP1\Monte_Carlo_DGP1_T_200_rho_0.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel C: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 50; rho = 0.5|---------
load("Results DGP1\Monte_Carlo_DGP1_T_50_rho_50.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel D: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 100; rho = 0.5|---------
load("Results DGP1\Monte_Carlo_DGP1_T_100_rho_50.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel E: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 200; rho = 0.5|---------
load("Results DGP1\Monte_Carlo_DGP1_T_200_rho_50.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel F: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 50; rho = 0.95|---------
load("Results DGP1\Monte_Carlo_DGP1_T_50_rho_95.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel G: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 100; rho = 0.95|---------
load("Results DGP1\Monte_Carlo_DGP1_T_100_rho_95.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel H: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 200; rho = 0.95|---------
load("Results DGP1\Monte_Carlo_DGP1_T_200_rho_95.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
Ta = array2table(MSFE_relative_rounded,'RowNames', row_headers,'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel I: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);


%% ========================================== |Results Table 2 Appendices (DGP1)| =================================================
clear;
disp('=============================================================== ');
disp('Appendices Table 3 (DGP2) ');
disp('=============================================================== ');
% ------------ |T = 50; rho = 0.00|---------
load("Results DGP2\Monte_Carlo_DGP2_T_50_rho_0.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel A: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 100; rho = 0.00|---------
load("Results DGP2\Monte_Carlo_DGP2_T_100_rho_0.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel B: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 200; rho = 0.00|---------
load("Results DGP2\Monte_Carlo_DGP2_T_200_rho_0.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel C: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 50; rho = 0.50|---------
load("Results DGP2\Monte_Carlo_DGP2_T_50_rho_50.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel D: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 100; rho = 0.50|---------
load("Results DGP2\Monte_Carlo_DGP2_T_100_rho_50.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel E: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 200; rho = 0.50|---------
load("Results DGP2\Monte_Carlo_DGP2_T_200_rho_50.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel F: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 50; rho = 0.95|---------
load("Results DGP2\Monte_Carlo_DGP2_T_50_rho_95.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel G: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 100; rho = 0.95|---------
load("Results DGP2\Monte_Carlo_DGP2_T_100_rho_95.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel H: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

% ------------ |T = 200; rho = 0.95|---------
load("Results DGP2\Monte_Carlo_DGP2_T_200_rho_95.mat");
MSFE_relative_rounded = round(MSFE_relative, 3);
% Define row and column headers
row_headers = {'TVP', 'Agnostic AVP (m=20)', 'Agnostic AVP (m=40)', ...
               'Agnostic AVP (m=60)', 'Targeted AVP (m=20)', ...
               'Targeted AVP (m=40)', 'Targeted AVP (m=60)'};
col_headers = {'beta_1', 'beta_2', 'beta_3', 'beta_4'};
Ta = array2table(MSFE_relative_rounded, 'RowNames', row_headers, 'VariableNames', col_headers);
disp('=============================================================== ');
disp(sprintf('Panel I: T = %d and rho = %.2f:', T, rho));
disp('=============================================================== ');
disp(Ta);

