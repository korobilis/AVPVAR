function [y,x,Z1,Z2,Z3,Z4,Z5,Z6,theta_t,beta_t,sigma_t] = outtvpsv_reg_dgp1(T,p,param)



%% ===========================| Structural Break DGP |===============================================



% **************************************************************************************************************************
% OUTTVP_REG_DGP  Code that generates data from a univariate time-varying parameter (tvp) regression with
% outliers. 
% The model is of the following form
%
%                   y[t]  =  x[t] beta[t] + e[t]
%               beta[t]  =  mu + 0.99I x (beta[t-1] - mu) + outlier + u[t]
%
% where e[t] ~ N(0, sigma), u[t] ~ N(0,(1/T)^2 x I), mu is the unconditional mean of the AR process for beta[t] 
% and "outlier" is a disrupting process
% **************************************************************************************************************************
% INPUTS:
%  T        Time series observations
%  p        Number of predictors   
%  param    Structure array that specified all other tuning parameters
%      param.rho     Correlation coefficient for predictors x
%      param.q       Number of predictors that have non-zero coefficients (i.e. they are important/"significant" for y
%      param.s       Generate the (Txp) variable s, which indexes which coefficients are zero (or not) in what periods
%      param.sigma0  Initial condition for the process for log(sigma[t])
%      param.theta0  Initial condition for the process for theta[t]
%      param.mu      Unconditional mean for the process for theta[t]
%
% OUTPUTS:
%  y        Generated time series following the sparse tvp regression
%  x        Generated right-hand side predictor variables
%  beta_t   Generated coefficients beta[t]
%  sigma_t  Regression variance
%
% **************************************************************************************************************************
% Written by Dimitris Korobilis on 11/06/2025
% University of Glasgow
% **************************************************************************************************************************

%% Check for INPUT arguments
if nargin == 0
    T   = 200;          % Time series observations
    p   = 8;            % Number of predictors
    rho = 0;            % Correlation between predictors is rho^|i-j| 
    
    V0 = log(.2);  % Regression variance
    theta0 = 2*(4.*rand(1,p)+1);  % Initial regression coefficients
    mu = theta0./2;     % Long-run (unconditional) mean of regression coefficients
else
    rho = param.rho;
    V0 = log(param.sigma0); 
    theta0 = param.theta0; 
    mu = param.mu;
end

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
    else
        theta_t(t,:) = mu +  0.999*(theta_t(t-1,:) - mu) + (1/(T^(2/4)))*randn(1,p);
    end
end
theta_t = theta_t(101:end,:);

% 3. Generate outlier process for the coefficients: Specific for 4 coeff
It = zeros(T,1); 
h1 = [round(6*T/8):round(8*T/8)];           % Outlier observation
w = 2 + 0.2*rand(length(h1),1);        % Outlier strength
It(h1,1) = w;
h2 = [round(5*T/8):round(8*T/8)];           % Outlier observation
w = 2 + 0.2*rand(length(h2),1);        % Outlier strength
It(h2,2) = w;
h3 = [round(4*T/8):round(8*T/8)];           % Outlier observation
w = 2 + 0.2*rand(length(h3),1);        % Outlier strength
It(h3,3) = w;
h4 = [round(3*T/8):round(8*T/8)];           % Outlier observation
w = 2 + 0.2*rand(length(h4),1);        % Outlier strength
It(h4,4) = w;

% Define total regression coefficients beta_t
beta_t = theta_t + It;



% Generate T x N matrix of standard normal random variables

% Non-informative Drivers
Z1 = [ones(T,1)  randn(T,20)  ];
Z2 = [ones(T,1)  randn(T,40)  ];
Z3 = [ones(T,1)  randn(T,60)  ];

% Informative Drivers
Z4 = [ones(T,1) It normalize(randn(T,10)) randn(T,10).*It(:,1) randn(T,10).*It(:,2) randn(T,10).*It(:,3) randn(T,10).*It(:,4) ];
Z5 = [ones(T,1) It normalize(randn(T,20)) randn(T,20).*It(:,1) randn(T,20).*It(:,2) randn(T,20).*It(:,3) randn(T,20).*It(:,4) ];
Z6 = [ones(T,1) It normalize(randn(T,30)) randn(T,30).*It(:,1) randn(T,30).*It(:,2) randn(T,30).*It(:,3) randn(T,30).*It(:,4) ];

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

