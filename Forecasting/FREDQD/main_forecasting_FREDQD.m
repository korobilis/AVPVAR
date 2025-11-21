%% Main Forecasting FRED-QD
clear; close all; clc;
% Add paths
addpath('functions');
addpath('data');
addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/
%===========================| USER INPUT |=================================
quant           = [10,25,50,75,90]./100;   
q10             = find(quant==0.1);  qmed = find(quant==0.5);  q90 = find(quant==0.9);
nq              = length(quant);   

%% LOOP OVER DIFFERENT LAG SPECIFICATIONS
p_values = [1 2 3 4];  % Different lag specifications to test

for p_idx = 1:length(p_values)
    
    p = p_values(p_idx);  % Current lag specification
    
    fprintf('\n========================================\n')
    fprintf('RUNNING FORECASTING ROUTINE FOR p = %d\n', p)
    fprintf('========================================\n\n')
    
    %% Please choose VAR specification: 
    select          = {'GDPC1', 'PCECTPI', 'FEDFUNDS'}; 
    stationY        = [5,5,2];
    station         = 0;  % 0: All instruments transformed as in MacCracken. 
    standarize      = 1;
    csv_in='current.csv';
    [data,dates,series,tcode] = fred_qd2(csv_in,select,station,stationY);
    data  = data(13:end,:); %missing values in instruments
    dates = dates(13:end,:); 
    datesCCMM = datenum(dates);  
    if standarize==1
        data = normalize(data);
    end

    %% Please select Instruments
    selectz          = {'S&P 500','S&P PE ratio','PERMITS','PPIIDC','BAA10YM','TOTALSLx','EXJPUSx','VIXCLSx','UMCSENTx'};
    [~, idx_selected] = ismember(selectz, series);
    Z = data(:, idx_selected);

    %% VAR inputs (p is now set in the loop)
    r               = 1;                     % Number of factors
    h               = 8;                     % Number of prediction steps ahead

    %% MCMC setting (total iterations are nsave+nburn, stored iterations are nsave/nthin)
    nsave    = 5000;           % Store nsave MCMC iterations from posteriors
    nburn    = 1000;           % Discard nburn iterations for convergence
    nthin    = 5;            % Save every nthin iteration to reduce correlation in the chain
    %% ===========================| LOAD DATA |=================================

    %% ======|Load series to extract factors from
    varnum = size(select,2);
    Y = data(:,1:varnum);
    T = size(Y,1);
    %% ========== Define the endogenous factor variable
    [FY]      = extract(zscore(Z),1);
    FY        = FY/chol(cov(FY)) - mean(FY/chol(cov(FY)));
    YFact     = [Y FY];

    %% =========== Start with  first_per% of observations for first period forecast
    first_per  = 0.5; 
    T_full     = T-h;
    T_thres    = round(first_per*T_full); 
    anumber    = T_full-T_thres+1;

    %Generate Matrices for storing results
    CQs        = zeros(anumber,nq,varnum,h,9);  % Save raw conditional quantiles   
    QScore10   = zeros(anumber,varnum,h,9);     % Quantile score 10%
    QScore90   = zeros(anumber,varnum,h,9);     % Quantile score 90%
    PL         = zeros(anumber,varnum,h,9);     % Predictive likelihood (based on kernel smoothing)
    PIT        = zeros(anumber,varnum,h,9);     % PITS (based on kernel smoothing)
    PE         = zeros(anumber,varnum,h,9);     % Prediction Error   
    MSPE       = zeros(anumber,varnum,h,9);     % Mean square prediction error (using median)
    MSFE       = zeros(anumber,varnum,h,9);     % Mean square forecast error (using mean)
    MAE        = zeros(anumber,varnum,h,9);     % Mean absolute error
    elapsedTime = zeros(anumber,9);              % Computation time
    invalid_runs_total = zeros(anumber,9);       % Invalid runs counter



    %% ======| Begin forecasting iterations
    for sample = T_thres:T_full %T_thres:T_full   sample = T_full
        disp('  ')           
        fprintf(['<strong>YOU ARE RUNNING SAMPLE ' num2str(sample) ' OF ' num2str(T_full) ' (p=' num2str(p) ')</strong> \n'] )    
        disp('  ')        

        yf_true = Y(sample+1:sample+h,:);  
        [FY]       = extract(zscore(Z(1:sample,:)),1);
         FY        = FY/chol(cov(FY)) - mean(FY/chol(cov(FY)));
        YFact      = [Y(1:sample,:) FY];

        %% =================== Model 1: AP-VAR
        tStart = tic;
        [yfore_save,beta_save,invalid_runs] = APVAR(Y(1:sample,:),Z(1:sample,1:end),h,p,r,nsave,nburn,nthin,'SV'); % SV: Stochastic Volatility. CV: Constant Volatility
        elapsedTime(sample-T_thres+1,1) = toc(tStart);
        invalid_runs_total(sample-T_thres+1,1) = invalid_runs;
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,1)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,1)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,1)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));         
                MAE(sample-T_thres+1,ivar,nfore,1)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,1)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,1)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end

     
        %% ========================= Model 2: CP-VAR
        tStart = tic;
        [yfore_save] = TVP_VAR(Y(1:sample,:),h,p,r,nsave,nburn,nthin,'CP','CL','CV'); 
        elapsedTime(sample-T_thres+1,2) = toc(tStart);
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,2)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,2)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,2)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,2)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,2)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,2)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end

         %% ============== Model 3: TVP-RW-EB
        tStart = tic;
        [yfore_save]   =  TVP_RW_EB(Y(1:sample,:),p,nsave,nburn,h);
        elapsedTime(sample-T_thres+1,3) = toc(tStart);
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,3)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,3)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,3)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));           
                MAE(sample-T_thres+1,ivar,nfore,3)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,3)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,3)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end

     %% ====================================== Model 4: CP-VAR SV
        tStart = tic;
        [yfore_save] = TVP_VAR(Y(1:sample,:),h,p,r,nsave,nburn,nthin,'CP','TVL-RW','SV');
        elapsedTime(sample-T_thres+1,4) = toc(tStart);
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,4)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,4)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,4)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,4)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,4)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,4)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end
     % 
     %% =================== Model 5: OLS
        tStart = tic;
        [yfore_save]   = BVAR_OLS_iter(Y(1:sample,:),p,h,nsave);
        elapsedTime(sample-T_thres+1,5) = toc(tStart);
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                          = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,5)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,5)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,5)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,5)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,5)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,5)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end


      %% ============== Model 6: SVO-t
        tStart = tic;
        [yfore_save]   = CCMM_SVO(Y(1:sample,:),sample,datesCCMM(1:sample),p,yf_true',h);  
        elapsedTime(sample-T_thres+1,6) = toc(tStart);

        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,6)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,6)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,6)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,6)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,6)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,6)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end
        
      %% =================== Model 7: FAVAR
        tStart = tic;
        [yfore_save]   = BVAR_OLS_iter(YFact(1:sample,:),p,h,nsave);
        elapsedTime(sample-T_thres+1,7) = toc(tStart);
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                          = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,7)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,7)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,7)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,7)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,7)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,7)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end

      %% ====================================== Model 8: FAVAR-SV
        tStart = tic;
        [yfore_save] = TVP_VAR(YFact(1:sample,:),h,p,r,nsave,nburn,nthin,'CP','TVL-RW','SV');
        elapsedTime(sample-T_thres+1,8) = toc(tStart);
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,8)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,8)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,8)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,8)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,8)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,8)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end

        %% ========================= Model 9: TVP_VAR
        tStart = tic;
        [yfore_save]  = TVP_VAR_FB(Y(1:sample,:),h,p,r,nsave,nburn,nthin);
        elapsedTime(sample-T_thres+1,9) = toc(tStart);
        invalid_runs_total(sample-T_thres+1,9) = invalid_runs;
        for ivar = 1:varnum
            for nfore = 1:h
                post_mean_quants                         = squeeze(quantile(yfore_save(ivar,nfore,:),quant));  
                CQs(sample-T_thres+1,:,ivar,nfore,9)     = (post_mean_quants);
                QScore10(sample-T_thres+1,ivar,nfore,9)  = (yf_true(nfore,ivar) - post_mean_quants(q10))*(double(yf_true(nfore,ivar)<=post_mean_quants(q10)) - quant(q10));
                QScore90(sample-T_thres+1,ivar,nfore,9)  = (yf_true(nfore,ivar) - post_mean_quants(q90))*(double(yf_true(nfore,ivar)<=post_mean_quants(q90)) - quant(q90));       
                MAE(sample-T_thres+1,ivar,nfore,9)       = abs(yf_true(nfore,ivar)-post_mean_quants(qmed));
                MSPE(sample-T_thres+1,ivar,nfore,9)      = (yf_true(nfore,ivar)-post_mean_quants(qmed))^2;
                MSFE(sample-T_thres+1,ivar,nfore,9)      = (yf_true(nfore,ivar)-squeeze(mean(yfore_save(ivar,nfore,:))))^2;
            end
        end

    end % End loop Forecasting

    % Save results for this lag specification
    filename = sprintf('oos_FREDQD_p%d.mat', p);
    save(filename);
    
    fprintf('\n========================================\n')
    fprintf('RESULTS SAVED TO: %s\n', filename)
    fprintf('========================================\n\n')

end % Loop over p values

fprintf('\n========================================\n')
fprintf('ALL LAG SPECIFICATIONS COMPLETED!\n')
fprintf('Files saved:\n')
for p_idx = 1:length(p_values)
    fprintf('  - oos_Final_p%d.mat\n', p_values(p_idx))
end
fprintf('========================================\n')