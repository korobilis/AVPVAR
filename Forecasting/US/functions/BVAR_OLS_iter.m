function [Y_forc,A_OLS,SIGMA_OLS,Residuals] = BVAR_OLS_iter(Y,p,h,ndraws)
% Based on Koop, Korobilis and Pettenuzzo (2015)
% =========================================================================
% This code estimates the Bayesian Compressed Vector Autoregression model
% (will expand)
% INPUTS         Y: TxM matrix of data,
%                p: number of lags,
%                h: forecast horizons,
%                ndraws: number of Monte Carlo draws
%                Note: ndraws = 1 gives point forecasts
% OUTPUT: Y_forc is the predictive density
% ======================| Set up and estimate VAR |========================
[T, M] = size(Y);       % Dimensions of VAR


% Take lags, and correct observations
Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
T = T-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
Y = Y(p+1:end,:);

% ==============| Start computation

% Storage matrices
Y_forc = zeros(M,h,ndraws);

X = [ones(T,1),Ylag];
    
% Analytical parameter posteriors
A_OLS = (X'*X)\(X'*Y);
SIGMA_OLS = (Y - X*A_OLS)'*(Y - X*A_OLS)./T;
Residuals = (Y - X*A_OLS);

        cSIGMA = chol(SIGMA_OLS);
        Y_pred = zeros(h,M);        
          

        for iter = 1:ndraws %  Very Inefficient
        X_fore = [Y(end,:) X(end,2:M*(p-1)+1)];  
        Y_pred(1,:) = [1 X_fore]*A_OLS   + randn(1,M)*cSIGMA;

        for ifore = 2:h
            X_fore = [Y_pred(ifore-1,:) X_fore(:,1:M*(p-1))];        
            Y_pred(ifore,:) = [1 X_fore]*A_OLS  +  randn(1,M)*cSIGMA;
        end

        Y_forc(:,:,iter)  =  Y_pred';
        end




