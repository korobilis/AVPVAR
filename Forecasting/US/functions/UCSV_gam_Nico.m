function [y_fitted, y_forecast, forecast_quantiles, tauhat,residuals, hhat, ghat, thetahat, thetastd, posterior_draws] = UCSV_gam_Nico(y, nloop, burnin,nthin, h_ahead)
% UCSV_GAM_NICO - Unobserved Components Stochastic Volatility Model with Forecasting
% ================================================================================
%  MODEL DESCRIPTION
% ================================================================================
% The Stock-Watson UCSV model decomposes y_t into:
% 1. MEASUREMENT EQUATION:
%    y_t = τ_t + ε_t,  where ε_t ~ N(0, exp(h_t))
%
%    - y_t: observed inflation at time t
%    - τ_t: unobserved trend/permanent component (this is the fitted value)
%    - ε_t: transitory component with stochastic volatility
%    - h_t: log-volatility of the transitory component
%
% 2. STATE EQUATION FOR TREND:
%    τ_t = τ_{t-1} + η_t,  where η_t ~ N(0, exp(g_t))
%    - g_t: log-volatility of the trend innovations
%
% 3. STOCHASTIC VOLATILITY EQUATIONS:
%    h_t = h_0 + ω_h * h̃_t
%    g_t = g_0 + ω_g * g̃_t
%
%    where h̃_t and g̃_t follow random walks:
%    h̃_t = h̃_{t-1} + ν_{h,t},  ν_{h,t} ~ N(0, 1)
%    g̃_t = g̃_{t-1} + ν_{g,t},  ν_{g,t} ~ N(0, 1)
%
% 4. FORECASTING:
%    For h-step ahead forecasts:
%    - τ_{T+h} = τ_T + Σ_{j=1}^h η_{T+j}
%    - h_{T+h} = h_0 + ω_h * (h̃_T + Σ_{j=1}^h ν_{h,T+j})
%    - g_{T+h} = g_0 + ω_g * (g̃_T + Σ_{j=1}^h ν_{g,T+j})
%    - y_{T+h} = τ_{T+h} + ε_{T+h}
%
% ================================================================================
% INPUTS:
%   y        - T×1 vector of observations (e.g., inflation rates)
%   nloop    - Number of MCMC iterations
%   burnin   - Number of burn-in iterations
%   h_ahead  - Forecast horizon (number of periods ahead)
%
% OUTPUTS:
%   y_fitted           - T×1 vector of fitted values (posterior mean of τ)
%   y_forecast         - h_ahead×1 vector of point forecasts (posterior mean)
%   forecast_quantiles - h_ahead×3 matrix of forecast quantiles [5%, 50%, 95%]
%   tauhat            - T×1 posterior mean of trend component
%   hhat              - T×1 posterior mean of exp(h_t/2)
%   ghat              - T×1 posterior mean of exp(g_t/2)
%   thetahat          - 4×1 posterior means [ω_h^2, ω_g^2, h_0, g_0]
%   thetastd          - 4×1 posterior std devs
%   posterior_draws   - Structure with all posterior draws
% ================================================================================

    %% ============================================================
    %% STEP 1: INITIALIZATION AND PRIOR SPECIFICATION
    %% ============================================================  
    % Extract data dimensions
    T = length(y);    
    % Prior hyperparameters 
    tau0 = 0;        % Prior mean for initial trend τ_0
    Vtau = 10;       % Prior variance for τ_0 (large = uninformative)
    
    b0 = 0;          % Prior mean for h_0
    Vh0 = 10;        % Prior variance for h_0
    Vh = 10;         % Prior variance for h̃_1
    Vomegah = 0.2;   % Prior variance for ω_h
    
    c0 = 0;          % Prior mean for g_0
    Vg0 = 10;        % Prior variance for g_0
    Vg = 10;         % Prior variance for g̃_1
    Vomegag = 0.2;   % Prior variance for ω_g
    
    %% ============================================================
    %% STEP 2: INITIALIZE MARKOV CHAIN
    %% ============================================================
    
    % Initialize volatility levels at sample variance
    sample_var = var(y);
    h0 = log(sample_var)/2;  % Initial h_0
    g0 = log(sample_var)/2;  % Initial g_0
    
    % Initialize volatility scaling parameters
    omegah2 = 0.2;  % Initial ω_h^2
    omegag2 = 0.2;  % Initial ω_g^2
    omegah = sqrt(omegah2);
    omegag = sqrt(omegag2);
    
    % Initialize standardized volatility processes
    htilde = zeros(T,1);  % h̃_t initialized at zero
    gtilde = zeros(T,1);  % g̃_t initialized at zero
    
    % Compute actual volatilities
    h = h0 + omegah*htilde;  % h_t = h_0 + ω_h * h̃_t
    g = g0 + omegag*gtilde;  % g_t = g_0 + ω_g * g̃_t
    
    % Create sparse difference matrix H for state equations
    H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    
    %% ============================================================
    %% STEP 3: STORAGE ALLOCATION
    %% ============================================================
    
    % Allocate storage for posterior draws
    n_saved = nloop;
    store_theta = zeros(n_saved, 4);     % [ω_h^2, ω_g^2, h_0, g_0]
    store_tau = zeros(n_saved, T);       % Trend component
    store_h = zeros(n_saved, T);         % Log-volatility of transitory
    store_g = zeros(n_saved, T);         % Log-volatility of trend
    store_forecast = zeros(n_saved, h_ahead);  % Forecasts
    
    % For posterior density estimation
    npts = 500;
    omh_grid = linspace(-1, 1, npts)';
    omg_grid = linspace(-1, 1, npts)';
    store_pomh = zeros(npts, 1);
    store_pomg = zeros(npts, 1);
    
    % Set random seeds for reproducibility
    rand('state', sum(100*clock));
    randn('state', sum(200*clock));
    
    %% ============================================================
    %% STEP 4: GIBBS SAMPLING MAIN LOOP
    %% ============================================================
    
    for loop = 1:n_saved
        
        %% --------------------------------------------------------
        %% 4.1: SAMPLE TREND COMPONENT τ
        %% --------------------------------------------------------
        % Conditional posterior: τ | y, h, g ~ N(τ̂, D_τ)
        % Precision matrix for state equation
        HinvStauH = H'*sparse(1:T,1:T,[1/Vtau*exp(-g(1)); exp(-g(2:end))])*H;      
        % Precision matrix for conditional posterior
        invDtau = HinvStauH + sparse(1:T,1:T,1./exp(h)); 
        % Prior mean contribution
        alptau = H\sparse(1,1,tau0,T,1);
        % Posterior mean
        tauhat_gibbs = invDtau\(HinvStauH*alptau + y./exp(h));
        % Sample from posterior
        tau = tauhat_gibbs + chol(invDtau,'lower')'\randn(T,1);
        
        %% --------------------------------------------------------
        %% 4.2: SAMPLE VOLATILITY (h̃, h_0, ω_h)
        %% --------------------------------------------------------        
        % Create transformed observations
        Ystar = log((y - tau).^2 + 0.0001);
        % Call specialized SV sampler
        [htilde, h0, omegah, omegahhat, Domegah] = ...
            SVRW_gam_Nico(Ystar, htilde, h0, omegah, b0, Vh0, Vh, Vomegah);
        % Update full volatility process
        h = h0 + omegah*htilde;
        omegah2 = omegah^2;
        
        %% --------------------------------------------------------
        %% 4.3: SAMPLE TREND VOLATILITY (g̃, g_0, ω_g)
        %% --------------------------------------------------------
        % Compute trend innovations
        trend_innovations = [(tau(1)-tau0)/sqrt(Vtau); tau(2:end)-tau(1:end-1)];       
        % Transform to linear model
        Ystar = log(trend_innovations.^2 + 0.0001);
        
        % Sample volatility parameters
        [gtilde, g0, omegag, omegaghat, Domegag] = SVRW_gam_Nico(Ystar, gtilde, g0, omegag, c0, Vg0, Vg, Vomegag);
        % Update full volatility process
        g = g0 + omegag*gtilde;
        omegag2 = omegag^2;
        
        %% --------------------------------------------------------
        %% 4.4: GENERATE h-STEP AHEAD FORECASTS
        %% --------------------------------------------------------
        if loop > burnin && mod(loop,nthin) == 0
            % Current values at time T
            tau_T = tau(T);
            htilde_T = htilde(T);
            gtilde_T = gtilde(T);
            
            % Storage for this iteration's forecasts
            y_forecast_iter = zeros(h_ahead, 1);
            
            for h_step = 1:h_ahead
                % Forecast volatilities h steps ahead
                % h̃_{T+h} = h̃_T + sum of h standard normal innovations
                htilde_Th = htilde_T + sqrt(h_step)*randn;
                h_Th = h0 + omegah*htilde_Th;
                
                % g̃_{T+h} = g̃_T + sum of h standard normal innovations
                gtilde_Th = gtilde_T + sqrt(h_step)*randn;
                g_Th = g0 + omegag*gtilde_Th;
                
                % Forecast trend: τ_{T+h} = τ_T + sum of h trend innovations
                % Each innovation η_t ~ N(0, exp(g_t))
                % For simplicity, we use the forecasted g_{T+h} for all innovations
                tau_Th = tau_T;
                for j = 1:h_step
                    % Evolve g for each step
                    gtilde_j = gtilde_T + sqrt(j)*randn;
                    g_j = g0 + omegag*gtilde_j;
                    tau_Th = tau_Th + sqrt(exp(g_j))*randn;
                end
                
                % Forecast observation: y_{T+h} = τ_{T+h} + ε_{T+h}
                % where ε_{T+h} ~ N(0, exp(h_{T+h}))
                y_forecast_iter(h_step) = tau_Th + sqrt(exp(h_Th))*randn;
            end
            
            % Store forecast for this iteration
            store_forecast(loop-burnin, :) = y_forecast_iter';
        end
        
        %% --------------------------------------------------------
        %% 4.5: PROGRESS REPORTING
        %% --------------------------------------------------------
        if mod(loop, 20000) == 0
            fprintf('Completed %d of %d iterations...\n', loop, nloop);
        end
        
        %% --------------------------------------------------------
        %% 4.6: STORE RESULTS AFTER BURN-IN
        %% --------------------------------------------------------
        if loop > burnin && mod(loop,nthin) == 0
            i = loop - burnin;
            
            % Store parameter draws
            store_tau(i,:) = tau';
            store_h(i,:) = h';
            store_g(i,:) = g';
            store_theta(i,:) = [omegah2, omegag2, h0, g0];
            
            % Accumulate posterior density estimates
            store_pomh = store_pomh + normpdf(omh_grid, omegahhat, sqrt(Domegah));
            store_pomg = store_pomg + normpdf(omg_grid, omegaghat, sqrt(Domegag));
        end
    end
    
    %% ============================================================
    %% STEP 5: POSTERIOR INFERENCE
    %% ============================================================
    
    % Compute posterior means and standard deviations
    thetahat = mean(store_theta)';  % [ω_h^2; ω_g^2; h_0; g_0]
    thetastd = std(store_theta)';
    
    % Posterior mean of trend (this is the fitted value)
    tauhat = mean(store_tau)';
    y_fitted = tauhat;  % Fitted values are the trend component
    
    % Posterior mean of volatilities (on standard deviation scale)
    hhat = mean(exp(store_h/2))';  % exp(h_t/2) = std dev of transitory
    ghat = mean(exp(store_g/2))';  % exp(g_t/2) = std dev of trend
    
    %% ============================================================
    %% STEP 6: COMPUTE FORECAST STATISTICS
    %% ============================================================
    % Point forecasts (posterior mean)
    y_forecast = mean(store_forecast)';
    % Forecast quantiles [5%, 50%, 95%]
    forecast_quantiles = quantile(store_forecast, [0.05, 0.50, 0.95])';
    
    %% ============================================================
    %% STEP 7: PREPARE OUTPUT
    %% ============================================================
    
    % Package posterior draws for output
    posterior_draws.tau = store_tau;
    posterior_draws.h = store_h;
    posterior_draws.g = store_g;
    posterior_draws.theta = store_theta;
    posterior_draws.pomh = store_pomh/(nloop-burnin);
    posterior_draws.pomg = store_pomg/(nloop-burnin);
    posterior_draws.forecasts = store_forecast;
    posterior_draws.tau_quantiles = quantile(store_tau, [0.05, 0.95])';
    posterior_draws.h_quantiles = quantile(store_h, [0.05, 0.95])';
    posterior_draws.g_quantiles = quantile(store_g, [0.05, 0.95])';
    
    %% ============================================================
    %% STEP 8: DISPLAY RESULTS
    %% ============================================================ 
    
    fprintf('\nModel Fit Statistics:\n');
    fprintf('--------------------\n');
    residuals = y - y_fitted;
    fprintf('RMSE:           %.4f\n', sqrt(mean(residuals.^2)));
    fprintf('MAE:            %.4f\n', mean(abs(residuals)));
    fprintf('R-squared:      %.4f\n', 1 - var(residuals)/var(y));
    
end