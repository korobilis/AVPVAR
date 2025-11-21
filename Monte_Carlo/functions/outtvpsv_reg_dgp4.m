function [y,x,Z1,Z2,Z3,Z4,Z5,Z6,theta_t,beta_t,sigma_t] = outtvpsv_reg_dgp4(T,p,param)

%% ==============| Contemporaneous structural break (Permanent) |==================

    rho = param.rho;
    V0 = log(param.sigma0); 
    theta0 = param.theta0; 
    mu = param.mu;


%% Generate sparse time-varying parameter model
% 1. First generate (possibly) correlated predictors
corr_x = zeros(p,p);          % Correlation matrix for predictors
for i = 1:p                   % do the lazy version, using for loops
    for j = 1:p           
        corr_x(i,j) = rho^(abs(i-j));
    end
end
x = zeros(T,p);
for t = 2:T+50
    x(t,:) = rho*x(t-1,:) + randn(1,p)*chol(corr_x);  % Generate RHS predictors
end
x = x(end-T+1:end,:);

% 2. Generate unrestricted regression coefficients theta_t
theta_t = zeros(T+100,p);
for t = 1:T+100
    if t == 1
        theta_t(t,:) = mu +  0.999*(theta0 - mu) + (1/(T^(2/4)))*randn(1,p);
    elseif t<round(0.75*T)
        theta_t(t,:) = mu +  0.999*(theta_t(t-1,:) - mu) + (1/(T^(2/4)))*randn(1,p);
    else
        theta_t(t,:) = mu +  0.999*(theta_t(t-1,:) - mu) + (10/(T^(2/4)))*randn(1,p);
    end
end
theta_t = theta_t(101:end,:);

% Define total regression coefficients beta_t
beta_t = theta_t;

% Generate T x N matrix of standard normal random variables
% Non-informative Drivers
Z1 = [ones(T,1)  randn(T,20)  ];
Z2 = [ones(T,1)  randn(T,40)  ];
Z3 = [ones(T,1)  randn(T,60)  ];

% Informative Drivers
breakpoint = round(0.75*T);
S = diag([ones(breakpoint,1) 10*ones(T-breakpoint+1,1)]);
Z4 = [ones(T,1) S*randn(T,20)];
Z5 = [ones(T,1) S*randn(T,40)];
Z6 = [ones(T,1) S*randn(T,60)];

% 4. Generate unrestricted regression coefficients sigma_t
V_t = zeros(T+100,1);
for t = 1:T+100
    if t == 1
        V_t(t,:) = V0 +  0.99*(V0 - V0) + (1/(T^(2/4)))*randn;
    else
        V_t(t,:) = V0 +  0.99*(V_t(t-1,:) - V0) + (1/(T^(2/4)))*randn;
    end
end
sigma_t = exp(V_t(101:end,:));

% 5. Generate dependent variable y
y = sum(x.*beta_t,2) + sqrt(sigma_t).*randn(T,1);

