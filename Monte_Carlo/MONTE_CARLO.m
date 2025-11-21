clear all;close all;clc;
addpath('functions')
addpath('MODELS')

%% -------------------------------PRELIMINARIES--------------------------------------
nMC    = 1000;                 % Number of DGPs to simulate
DGP    = 'DGP1';               % Choose DGP: 'DGP1' 'DGP2'
ngibbs = 10000;                % Number of Gibbs sampler iterations
nburn  = 0.1*ngibbs;           % Number of iterations to discard
nthin  = 10;                   % Thin to reduce correlation in Gibbs chain
                      
% Define parameter grids
T_values = [50 100 200];       % Number of time series observations
rho_values = [0, 0.5, 0.95];   % Correlation among predictors
p = 4;                         % Number of predictors

%% ----------------------------- START PARAMETER LOOPS --------------------------------

% Initialize results storage structure
results = struct();

for i_T = 1:length(T_values)
    T = T_values(i_T);
    
    for i_rho = 1:length(rho_values)
        rho = rho_values(i_rho);
        
        fprintf('=================================================================\n');
        fprintf('Starting simulation for DGP=%s, T=%d, rho=%.2f\n', DGP, T, rho);
        fprintf('=================================================================\n');
        
        tic;
        
        %% ----------------------------------PREP DATA----------------------------------------   
        theta0 = 4.*rand(1,p)-2;  mu = theta0;
        param.theta0=theta0; param.sigma0 = 0.2; param.rho = rho; param.mu = mu; param.smoothness = 0.1; 
        
        %% ==============| Start Monte Carlo |================================
        BETA   = zeros(T,p,nMC,8);    
        SIGMA  = zeros(T,nMC,8);
        SFE    = zeros(nMC,p,7);      
        
        for iMC=1:nMC
            tic_mc = tic;
            
            %% Generate artificial data:
            switch DGP

                case 'DGP1'
                   [y,X,Z1,Z2,Z3,Z4,Z5,Z6,theta_sim,beta_sim,sigma_sim] = outtvpsv_reg_dgp5(T,p,param);  
                case 'DGP2'
                   [y,X,Z1,Z2,Z3,Z4,Z5,Z6,theta_sim,beta_sim,sigma_sim] = nonlinear_tvp_dgp(T);   
            end

            %% Estimate models
            % 1) KR1 algorithm (TVP benchmark)
            [beta_save1,sigma_save1] = KR1(y,X,ngibbs,nburn,3);
            
            % Agnostic Drivers
            % 2.1) APVAR
            [beta_save2,sigma_save2] = AP_RW(y,X,Z1,ngibbs,nburn);

            % 2.2) APVAR
            [beta_save3,sigma_save3] = AP_RW(y,X,Z2,ngibbs,nburn);

            % 2.3) APVAR
            [beta_save4,sigma_save4] = AP_RW(y,X,Z3,ngibbs,nburn);

            % Targeted Drivers
            % 3.1) APVAR
            [beta_save5,sigma_save5] = AP_RW(y,X,Z4,ngibbs,nburn);

            % 3.2) APVAR
            [beta_save6,sigma_save6] = AP_RW(y,X,Z5,ngibbs,nburn);

            % 3.3) APVAR
            [beta_save7,sigma_save7] = AP_RW(y,X,Z6,ngibbs,nburn);

            %% Do thinning
            beta_save1 = beta_save1(1:nthin:end,:,:);   sigma_save1 = sigma_save1(1:nthin:end,:);
            beta_save2 = beta_save2(1:nthin:end,:,:);   sigma_save2 = sigma_save2(1:nthin:end,:);
            beta_save3 = beta_save3(1:nthin:end,:,:);   sigma_save3 = sigma_save3(1:nthin:end,:);
            beta_save4 = beta_save4(1:nthin:end,:,:);   sigma_save4 = sigma_save4(1:nthin:end,:);
            beta_save5 = beta_save5(1:nthin:end,:,:);   sigma_save5 = sigma_save5(1:nthin:end,:);
            beta_save6 = beta_save6(1:nthin:end,:,:);   sigma_save6 = sigma_save6(1:nthin:end,:);
            beta_save7 = beta_save7(1:nthin:end,:,:);   sigma_save7 = sigma_save7(1:nthin:end,:);
            
            %% Save all coefficients
            BETA(:,:,iMC,1) = beta_sim;           SIGMA(:,iMC,1) = sigma_sim;
            BETA(:,:,iMC,2) = mean(beta_save1);   SIGMA(:,iMC,2) = mean(sigma_save1);
            BETA(:,:,iMC,3) = mean(beta_save2);   SIGMA(:,iMC,3) = mean(sigma_save2);
            BETA(:,:,iMC,4) = mean(beta_save3);   SIGMA(:,iMC,4) = mean(sigma_save3);
            BETA(:,:,iMC,5) = mean(beta_save4);   SIGMA(:,iMC,5) = mean(sigma_save4);
            BETA(:,:,iMC,6) = mean(beta_save5);   SIGMA(:,iMC,6) = mean(sigma_save5);
            BETA(:,:,iMC,7) = mean(beta_save6);   SIGMA(:,iMC,7) = mean(sigma_save6); 
            BETA(:,:,iMC,8) = mean(beta_save7);   SIGMA(:,iMC,8) = mean(sigma_save7);  

            SFE(iMC,:,1) = mean((beta_sim - squeeze(mean(beta_save1))).^2);
            SFE(iMC,:,2) = mean((beta_sim - squeeze(mean(beta_save2))).^2);
            SFE(iMC,:,3) = mean((beta_sim - squeeze(mean(beta_save3))).^2);
            SFE(iMC,:,4) = mean((beta_sim - squeeze(mean(beta_save4))).^2);
            SFE(iMC,:,5) = mean((beta_sim - squeeze(mean(beta_save5))).^2);
            SFE(iMC,:,6) = mean((beta_sim - squeeze(mean(beta_save6))).^2);
            SFE(iMC,:,7) = mean((beta_sim - squeeze(mean(beta_save7))).^2);

            if mod(iMC, 5) == 0
                fprintf('DGP=%s, T=%d, rho=%.2f: Completed %d out of %d MC iterations\n', DGP, T, rho, iMC, nMC);
            end
            toc(tic_mc);
        end

        %% Calculate results for this configuration
        models = {'KR1';'Agnostic AVP-VAR (m=small)';'Agnostic AVP-VAR (m=medium)';'Agnostic AVP-VAR (m=large)';'Targeted AVP-VAR (m=small)';'Targeted AVP-VAR (m=medium)';'Targeted AVP-VAR (m=large)'};

        MSFE = [mean(SFE(:,:,1),1); mean(SFE(:,:,2),1); mean(SFE(:,:,3),1); mean(SFE(:,:,4),1); ...
                mean(SFE(:,:,5),1); mean(SFE(:,:,6),1); mean(SFE(:,:,7),1)];
        MSFE_relative = MSFE ./ mean(SFE(:,:,1),1);

        config_name = sprintf('%s_T_%d_rho_%d', DGP, T, round(rho*100));
        
        results.(config_name).T = T;
        results.(config_name).rho = rho;
        results.(config_name).BETA = BETA;
        results.(config_name).SIGMA = SIGMA;
        results.(config_name).SFE = SFE;
        results.(config_name).MSFE = MSFE;
        results.(config_name).MSFE_relative = MSFE_relative;
        results.(config_name).models = models;
        results.(config_name).nMC = nMC;
        results.(config_name).DGP = DGP;
        results.(config_name).param = param;
        
        %% Save individual configuration results
        filename = sprintf('Monte_Carlo_%s_T_%d_rho_%d.mat', DGP, T, round(rho*100));
        
        % Save individual file for this configuration
        save(filename, 'BETA', 'SIGMA', 'SFE', 'MSFE', 'MSFE_relative', 'models', ...
             'T', 'rho', 'nMC', 'DGP', 'param', 'ngibbs', 'nburn', 'nthin');
        
        
        fprintf('Results saved for T=%d, rho=%.2f in file: %s\n', T, rho, filename);
        
        elapsed_time = toc;
        fprintf('Completed T=%d, rho=%.2f in %.2f seconds\n', T, rho, elapsed_time);
        fprintf('-----------------------------------------------------------------\n');

         save(sprintf('Monte_Carlo_%s_All_Configurations.mat', DGP), 'results', 'T_values', 'rho_values', ...
             'nMC', 'DGP', 'ngibbs', 'nburn', 'nthin');
        
    end 
end 



