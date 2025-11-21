%% INSAMPLE models
clear; close all; clc;
%% ======| Add paths
addpath('data')
addpath('functions')
%% ===========================| USER INPUT |=================================
quant = [10,25,50,75,90]./100;   
q10 = find(quant==0.1);  qmed = find(quant==0.5);  q90 = find(quant==0.9);
nq = length(quant);   

%% ======| Please choose VAR specification:
select     = {'INDPRO', 'PCEPI', 'FEDFUNDS'};  
station    = 2;            
tcode      = (station==1)*[ 5, 5, 1] + (station==2)*[ 5, 5, 2];
standard   = 1;      

%% ======|Please choose Instruments Z:
selectZ    = {'USEPU','GPR','GECON','GZ','EBP','JLN12','JLNF12','UMCSENT','ISENT','GFC','ADS','EMV'}; %1985 and after Instrumets
stationZ   = 1;            
tcodeZ     = (stationZ==1)*[1,1,1,1,1,1,1,1,1,1,1,1] + (stationZ==2)*[2,2,2,2,2,2,2,2,2,2,2,2];
standardZ  = 1;  

%% ======| VAR inputs
p          = 1;                     % Number of lags
r          = 1;                     % Number of factors
h          = 1;                    % Number of prediction steps ahead

%% ======| MCMC setting (total iterations are nsave+nburn, stored iterations are nsave/nthin)
nsave    = 95000;         % Store nsave MCMC iterations from posteriors 95000
nburn    = 5000;         % Discard nburn iterations for convergence 5000
nthin    = 5;            % Save every nthin iteration to reduce correlation in the chain
%% ===========================| LOAD DATA |=================================
%% ======|Load series to extract factors from
[Y,Z,T,varnum,data,dataZ,dates,varnames] = load_data_Instruments(select,selectZ,tcode,tcodeZ,standard,standardZ);
 
%% ======| Start with  first_per% of observations for first period forecast
T_full    = T-h;
sample = T_full   
model_names={'Model 1: AP-VAR','Model 2: TVP-VAR-EB','Model 3: TVP-VAR-FB'};

     
%% =================== Model 1: AP-VAR
[yfore_save,beta_save1,~,~,residuals_save1] = APVAR(Y(1:sample,:),Z(1:sample,1:end),h,p,r,nsave,nburn,nthin,'SV');
%% ============== Model 2:  TVP-VAR-EB
[residuals_save2,~,beta_save2] = TVP_RW_EB_residuals(Y(1:sample,:),p,nsave,nburn,nthin);
%% ============== Model 3:  TVP-VAR-FB
[yfore_save3,beta_save3,invalid_runs,~,residuals_save3] = TVP_VAR_FB(Y(1:sample,:),h,p,r,nsave,nburn,nthin);
%% ============== UCSV
[y_fitted(:,1), y_forecast(:,1),~, tauhat(:,1),residualsUCSV(:,1)] = UCSV_gam_Nico(Y(1:sample,1), nsave, nburn,nthin, h);
[y_fitted(:,2), y_forecast(:,2),~, tauhat(:,2),residualsUCSV(:,2)] = UCSV_gam_Nico(Y(1:sample,2), nsave, nburn,nthin, h);
[y_fitted(:,3), y_forecast(:,3),~, tauhat(:,3),residualsUCSV(:,3)] = UCSV_gam_Nico(Y(1:sample,3), nsave, nburn,nthin, h);

% Save
%save('insample_p1.mat')





