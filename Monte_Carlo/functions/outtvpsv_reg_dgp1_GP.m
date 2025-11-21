function [y,x,Z1,Z2,Z3,Z4,Z5,Z6,theta_t,beta_t,sigma_t] = outtvpsv_reg_dgp1_GP(T,p,param)
%% 
% This DGP does GP parameters, and also includes a structural break
%%



V0 = log(param.sigma0); 

%% 1. Generate predictors 
corr_x = zeros(p,p);
for i = 1:p
    for j = 1:p           
        corr_x(i,j) = param.rho^(abs(i-j));
    end
end
x = zeros(T,p);
for t = 2:T+50
    x(t,:) = param.rho*x(t-1,:) + randn(1,p)*chol(corr_x);
end
x = x(end-T+1:end,:);

%% 2. Generate time-varying parameters using different non-parametric approaches
theta_t = generate_gp_tvp(T, p, param);
% 3. Generate outlier process for the coefficients

h = [round(T/2):round(4*T/4)];           % BREAK
w = 2 + 0.2*rand(length(h),1);       % Outlier strength
It = zeros(T,1); It(h,:) = w;
% Define total regression coefficients beta_t
beta_t = theta_t + It;



%% 4. Generate auxiliary variables Z1, Z2
% Non-informative Drivers
Z1 = [ones(T,1)  randn(T,20)  ];
Z2 = [ones(T,1)  randn(T,40)  ];
Z3 = [ones(T,1)  randn(T,60)  ];

% Informative Drivers
Z4 = [ones(T,1) It normalize(randn(T,10)) randn(T,10).*It];
Z5 = [ones(T,1) It normalize(randn(T,20)) randn(T,20).*It];
Z6 = [ones(T,1) It normalize(randn(T,30)) randn(T,30).*It];






%% 5. Generate unrestricted regression coefficients sigma_t
V_t = zeros(T+100,1);
for t = 1:T+100
    if t == 1
        V_t(t,:) = V0 +  0.99*(V0 - V0) + (1/(T^(2/4)))*randn;
    else
        V_t(t,:) = V0 +  0.99*(V_t(t-1,:) - V0) + (1/(T^(2/4)))*randn;
    end
end
sigma_t = exp(V_t(end-T+1:end,:));

%% 6. Generate dependent variable
y = sum(x.*beta_t,2) + sqrt(sigma_t).*randn(T,1);

end

%% ==================== GP
% Gaussian Process with RBF kernel
function theta_t = generate_gp_tvp(T, p, param)

    time_grid = (1:T)'/T;
    theta_t = zeros(T, p);
    
    % GP hyperparameters
    length_scale = param.l; % l
    sigma_f = 1;  % Signal variance
    sigma_n = 0.001; % Noise variance
    
    % Compute covariance matrix
    K = zeros(T, T);
    for i = 1:T
        for j = 1:T
            K(i,j) = sigma_f^2 * exp(-0.5 * ((time_grid(i) - time_grid(j))/length_scale)^2);
        end
    end
    K = K + sigma_n^2 * eye(T);  % Add noise term
    
    % Cholesky decomposition for sampling
    L = chol(K, 'lower');
    
    for j = 1:p
        % Sample from GP
        z = randn(T, 1);
        f_sample = L * z;
        
        % Add mean function
        if isfield(param, 'theta0')
            theta_t(:, j) = param.theta0(j) + f_sample;
        else
            theta_t(:, j) = f_sample;
        end
    end
end

