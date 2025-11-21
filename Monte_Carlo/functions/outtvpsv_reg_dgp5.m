function [y,x,Z1,Z2,Z3,Z4,Z5,Z6,theta_t,beta_t,sigma_t] = outtvpsv_reg_dgp5(T,p,param)

%% ==============| Contemporaneous structural break (Jump)|==================

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
for t = 2:T+100
    x(t,:) = rho*x(t-1,:) + randn(1,p)*chol(corr_x);  % Generate RHS predictors
end
x = x(end-T+1:end,:);

% 2. Generate unrestricted regression coefficients theta_t
theta_t = zeros(T+100,p);
for t = 1:T+100
    if t == 1
        theta_t(t,:) = mu +  0.999*(theta0 - mu) + (1/(T^(2/4)))*randn(1,p);
    else
        theta_t(t,:) = mu +  0.999*(theta_t(t-1,:) - mu) + (1/(T^(2/4)))*randn(1,p);
    end
end
theta_t = theta_t(101:end,:);

% 3. Generate outlier process for the coefficients
h = [round(2*T/3):round(2*T/3)+5];          % Outlier observation h = [round(2*T/3):round(3*T/4)];   [floor(5*T/6):ceil(7*T/8)]; 
w = 6 + 0.2*rand(length(h),1);              % Outlier strength: w = 2 + 0.2*rand(length(h),1);   
It = zeros(T,1); It(h,:) = w;

% Define total regression coefficients beta_t
beta_t = theta_t + It;



% Generate T x N matrix of standard normal random variables

% Non-informative Drivers
Z1 = [ones(T,1)  randn(T,20)  ];
Z2 = [ones(T,1)  randn(T,40)  ];
Z3 = [ones(T,1)  randn(T,60)  ];



% Informative Drivers
Z4 = [ones(T,1) It normalize(randn(T,10)) randn(T,10).*It];
Z5 = [ones(T,1) It normalize(randn(T,20)) randn(T,20).*It];
Z6 = [ones(T,1) It normalize(randn(T,30)) randn(T,30).*It];



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

