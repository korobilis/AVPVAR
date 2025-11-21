function [y,x,Z1,Z2,Z3,Z4,Z5,Z6,theta_t,beta_t,sigma_t,regime_t] = nonlinear_tvp_dgp(T,p,param)

% Nonlinear TVP DGP with Regime Switching and Threshold Effects
% This DGP features:
% 1. Smooth transition between regimes based on lagged y
% 2. Time-varying coefficients that depend on regime and smooth transitions
% 3. Heteroskedastic errors that depend on regime
% 4. Informative instruments that capture regime information

%% Check for INPUT arguments
if nargin == 1
    %T = 200;
    p = 4;
    % Default parameters
    param.rho = 0.3;
    param.sigma0 = 0.2;
    param.theta0 = [1.5, 0.8, -0.4, 0.2];
    param.mu = param.theta0 ./ 2;
    param.gamma = 2.0;  % Transition smoothness parameter
    param.c = 0.0;      % Threshold parameter
else
    if ~isfield(param, 'gamma'), param.gamma = 2.0; end
    if ~isfield(param, 'c'), param.c = 0.0; end
end

rho = param.rho;
V0 = log(param.sigma0);
theta0 = param.theta0;
mu = param.mu;
gamma = param.gamma;
c = param.c;

%% 1. Generate correlated predictors (same as before)
corr_x = zeros(p,p);
for i = 1:p
    for j = 1:p           
        corr_x(i,j) = rho^(abs(i-j));
    end
end
x = zeros(T,p);
for t = 2:T+50
    x(t,:) = rho*x(t-1,:) + randn(1,p)*chol(corr_x);
end
x = x(end-T+1:end,:);

%% 2. Generate regime-switching nonlinear TVP process
theta_t = zeros(T,p);
regime_t = zeros(T,1);  % Regime indicator (0 or 1)
S_t = zeros(T,1);       % Smooth transition function
y = zeros(T,1);

% Initialize
y(1) = 0;
theta_t(1,:) = theta0;

for t = 2:T
    % Smooth transition function based on lagged y
    if t > 2
        S_t(t) = 1 / (1 + exp(-gamma * (y(t-1) - c)));
    else
        S_t(t) = 0.5;  % Initial value
    end
    
    % Regime indicator (for instruments)
    regime_t(t) = S_t(t) > 0.7;  % High regime when S_t > 0.7
    
    % Time-varying parameters with heterogeneous nonlinear evolution
    for j = 1:p
        % Base AR(1) evolution with parameter-specific persistence
        persistence = 0.95 + 0.04*(j-1)/(p-1);  % 0.95 to 0.99 across parameters
        theta_base = mu(j) + persistence * (theta_t(t-1,j) - mu(j));
        
        if t > 2
            % Parameter 1: Primarily regime-driven (strong regime effects)
            if j == 1
                regime_drift = 0.4 * S_t(t) * (2*regime_t(t) - 1);
                interaction_effect = 0.05 * tanh(x(t-1,j)) * S_t(t);
                stress_measure = abs(y(t-1)) + 0.5*abs(y(t-2));
                threshold_effect = 0.1 * (stress_measure > 1.5) * (2*regime_t(t) - 1);
            
            % Parameter 2: Interaction-driven (strong predictor interactions)
            elseif j == 2
                regime_drift = 0.15 * S_t(t) * (2*regime_t(t) - 1);
                interaction_effect = 0.25 * tanh(2*x(t-1,j)) * S_t(t);
                stress_measure = abs(y(t-1)) + 0.5*abs(y(t-2));
                threshold_effect = 0.05 * (stress_measure > 1.2) * (2*regime_t(t) - 1);
                
            % Parameter 3: Threshold-driven (strong discontinuous effects)
            elseif j == 3
                regime_drift = 0.1 * S_t(t) * (2*regime_t(t) - 1);
                interaction_effect = 0.08 * tanh(x(t-1,j)) * S_t(t);
                stress_measure = abs(y(t-1)) + 0.8*abs(y(t-2));
                threshold_effect = 0.35 * (stress_measure > 1.0) * (2*regime_t(t) - 1);
                
            % Parameter 4: Mixed effects (moderate all effects)
            else
                regime_drift = 0.2 * S_t(t) * (2*regime_t(t) - 1);
                interaction_effect = 0.12 * tanh(x(t-1,j)) * S_t(t);
                stress_measure = abs(y(t-1)) + 0.3*abs(y(t-2));
                threshold_effect = 0.15 * (stress_measure > 1.3) * (2*regime_t(t) - 1);
            end
            
            theta_t(t,j) = theta_base + regime_drift + interaction_effect + threshold_effect + ...
                          (0.03 + 0.04*j/p)/sqrt(T) * randn;  % Increasing variance by parameter
        else
            theta_t(t,j) = theta_base + (0.03 + 0.04*j/p)/sqrt(T) * randn;
        end
    end
    
    % Generate y with regime-dependent variance
    base_variance = exp(V0 + 0.99*(log(0.2) - V0) + (0.05/sqrt(T))*randn);
    regime_variance = base_variance * (1 + 0.5 * S_t(t));  % Higher variance in high regime
    
    y(t) = sum(x(t,:) .* theta_t(t,:)) + sqrt(regime_variance) * randn;
end

beta_t = theta_t;  % For compatibility
sigma_t = ones(T,1) * param.sigma0;  % Simplified for now

%% 3. Generate Instruments
% First create the fundamental signals that drive parameter changes
regime_signal = zeros(T,1);
interaction_signal = zeros(T,1);  
threshold_signal = zeros(T,1);

for t = 2:T
    if t > 3
        % Direct linear signals that correspond to parameter drivers
        regime_signal(t) = S_t(t) * (2*regime_t(t) - 1);  % This directly affects theta_t
        interaction_signal(t) = S_t(t) * mean(tanh(x(t-1,:)));  % Interaction component
        stress_measure = abs(y(t-1)) + 0.5*abs(y(t-2));
        threshold_signal(t) = (stress_measure > 1.0) * (2*regime_t(t) - 1);  % Threshold component
    end
end

% Non-informative Drivers (pure noise)
Z1 = [ones(T,1)  randn(T,20)];   
Z2 = [ones(T,1)  randn(T,40)];   
Z3 = [ones(T,1)  randn(T,60)];

% Informative Drivers - LINEAR combinations of the structural signals
% Following your outlier DGP pattern: [1, signal, noise, noise.*signal]

Z4 = [ones(T,1), regime_signal, randn(T,9), randn(T,10).*regime_signal];

Z5 = [ones(T,1), regime_signal, interaction_signal, threshold_signal, ...
      randn(T,17), randn(T,10).*regime_signal, randn(T,10).*interaction_signal];

Z6 = [ones(T,1), regime_signal, interaction_signal, threshold_signal, ...
      randn(T,27), randn(T,10).*regime_signal, randn(T,10).*interaction_signal, ...
      randn(T,10).*threshold_signal];

end