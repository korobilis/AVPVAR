function [residuals_save, fitted_save, Btdraw_save] = TVP_RW_EB_residuals(Y,p,nrep,nburn,nthin)

t=size(Y,1); % t is the time-series observations of Y
M=size(Y,2); % M is the dimensionality of Y

% Number of factors & lags:
tau = 40; % tau is the size of the training sample
numa = M*(M-1)/2; % Number of lower triangular elements of A_t (other than 0's and 1's)

%% ===================================| VAR EQUATION |==============================
% Generate lagged Y matrix. This will be part of the X matrix
ylag = mlag2(Y,p); % Y is [T x M]. ylag is [T x (Mp)]
% Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = ylag(p+tau+1:t,:);
K = M + p*(M^2); % K is the number of elements in the state vector beta_t

% Calculate coefficients per equation
K_per_eq = K / M;  % This should be: 1 + p*M (constant + p lags of M variables)

% Verify this makes sense
if mod(K, M) ~= 0
    error('Total coefficients K must be divisible by number of equations M');
end

% Create Z_t matrix (regressors)
Z = zeros((t-tau-p)*M,K);
for i = 1:t-tau-p
    ztemp = eye(M);
    for j = 1:p        
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];  
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end

% Redefine VAR variables y after the training sample
y = Y(tau+p+1:t,:)';
% Time series observations
t=size(y,2);   % t is now T - p - tau

% Calculate number of saved iterations with thinning
nsave = floor(nrep/nthin);

% MODIFIED: Initialize storage with 4D structure to match BIVAR
% [number of draws saved] × [T] × [coefficients per equation] × [number of equations]
residuals_save = zeros(nsave, t, M);  % [draws × time × equations]
fitted_save = zeros(nsave, t, M);     % [draws × time × equations]
Btdraw_save = zeros(nsave, t, K_per_eq, M);  % 4D structure

%----------------------------PRELIMINARIES---------------------------------

%========= PRIORS:
% To set up training sample prior a-la Primiceri, use the following subroutine
[B_OLS,VB_OLS,A_OLS,sigma_OLS,VA_OLS]= ts_prior(Y,tau,M,p);

% Set some hyperparameters here (see page 831 of Primiceri, end of section 4.1)
k_Q = 0.01; 
k_S = 0.1;
k_W = 0.01;

% We need the sizes of some matrices as prior hyperparameters (see page
% 831 again, lines 2-3 and line 6)
sizeW = M; % Size of matrix W
sizeS = 1:M; % Size of matrix S

%-------- Now set prior means and variances (_prmean / _prvar)
% These are the Kalman filter initial conditions for the time-varying
% parameters B(t), A(t) and (log) SIGMA(t). These are the mean VAR
% coefficients, the lower-triangular VAR covariances and the diagonal
% log-volatilities, respectively 
% B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;
% A_0 ~ N(A_OLS, 4Var(A_OLS)) starting values for kalman recursion for
% parameters L
A_0_prmean = A_OLS;
A_0_prvar = 4*VA_OLS;
% log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean = sigma_OLS;
sigma_prvar = eye(M);

% Set priors for Q, S, W
% Note that for IW distribution I keep the _prmean/_prvar notation....
% Q is the covariance of B(t), S is the covariance of A(t) and W is the
% covariance of (log) SIGMA(t)
% Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = ((k_Q)^2)*tau*VB_OLS;
Q_prvar = tau;
% W ~ IG(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = ((k_W)^2)*ones(M,1);
W_prvar = 2;
% S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
S_prmean = cell(M-1,1);
S_prvar = zeros(M-1,1);
ind = 1;
for ii = 2:M
    % S is block diagonal as in Primiceri (2005)
    S_prmean{ii-1} = ((k_S)^2)*(1 + sizeS(ii-1))*VA_OLS(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind);
    S_prvar(ii-1) = 1 + sizeS(ii-1);
    ind = ind + ii;
end

%========= INITIALIZE MATRICES:
% Specify covariance matrices for measurement and state equations
consQ = 0.0001;
consS = 0.0001;
consH = 0.01;
consW = 0.0001;
% this is matrix Sigma stucked one after another
Ht = kron(ones(t,1),consH*eye(M));   % Initialize Htdraw, a draw from the VAR covariance matrix
Htchol = kron(ones(t,1),sqrt(consH)*eye(M)); % Cholesky of Htdraw defined above
Qdraw = consQ*eye(K);   % Initialize Qdraw, a draw from the covariance matrix Q
Sdraw = consS*eye(numa);  % Initialize Sdraw, a draw from the covariance matrix S
Sblockdraw = cell(M-1,1); % ...and then get the blocks of this matrix (see Primiceri)
ijc = 1;
for jj=2:M
    Sblockdraw{jj-1} = Sdraw(((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc);
    ijc = ijc + jj;
end
Wdraw = consW*ones(M,1);    % Initialize Wdraw, a draw from the covariance matrix W
Btdraw = zeros(K,t);     % Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Atdraw = zeros(numa,t);  % Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
Sigtdraw = zeros(t,M);                % Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
sigt = kron(ones(t,1),0.01*eye(M));   % Matrix of the exponent of Sigtdraws (SIGMA(t))
statedraw = 5*ones(t,M);              % initialize the draw of the indicator variable 
                                      % (of 7-component mixture of Normals approximation)
Zs = kron(ones(t,1),eye(M));

%----------------------------- END OF PRELIMINARIES ---------------------------

%% ====================================== START SAMPLING ========================================
tic; fprintf('Primiceri  ')

disp('Number of iterations');
counter = 0;
for irep = 1:nrep + nburn    % GIBBS iterations starts here
    % Print iterations
    if mod(irep,1000) == 0
        disp(irep);toc;
    end

    %%   STEP I: Sample B from p(B|y,A,Sigma,V) (Drawing coefficient states, pp. 844-845)
    draw_beta
    
    %%   STEP II: Draw A(t) from p(At|y,B,Sigma,V) (Drawing coefficient states, p. 845)
    draw_alpha
    
    %%   STEP III: Draw diagonal VAR covariance matrix log-SIGMA(t)
    draw_sigma
    
    %% Create the VAR covariance matrix H(t). It holds that:
    %           A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '
    Ht = zeros(M*t,M);
    Htsd = zeros(M*t,M);
    for i = 1:t % for each t construct sigma
        inva = inv(capAt((i-1)*M+1:i*M,:));
        stem = diag(sigt(i,:));
        Hsd = inva*stem;
        Hdraw = Hsd*Hsd';
        Ht((i-1)*M+1:i*M,:) = Hdraw;  % H(t)
        Htsd((i-1)*M+1:i*M,:) = Hsd;  % Cholesky of H(t)
    end
    
    %% ----------------------------IN-SAMPLE RESIDUALS AFTER BURN-IN DRAWS  -----------------
    if irep > nburn && mod(irep-nburn, nthin) == 0
        counter = counter + 1;
        
        % Storage for current iteration
        Y_residuals = zeros(M, t);
        Y_fitted = zeros(M, t);
        
        % Loop through all time periods to compute residuals
        for t_idx = 1:t
            % Extract the corresponding regressors from Z matrix for time t_idx
            Z_t = Z((t_idx-1)*M+1:t_idx*M, :);  % M x K matrix for time t_idx
            B_t = Btdraw(:, t_idx);              % K x 1 coefficient vector for time t_idx
            
            % Compute fitted values: y_hat_t = Z_t * B_t
            y_fitted = Z_t * B_t;               % M x 1 fitted values
            Y_fitted(:, t_idx) = y_fitted;
            
            % Compute residuals: e_t = y_t - y_hat_t
            residual = y(:, t_idx) - y_fitted;  % M x 1 residuals
            Y_residuals(:, t_idx) = residual;
        end
        
        % MODIFIED: Reshape and store Btdraw in 4D format
        % Btdraw is currently [K x T], need to reshape to [T x K_per_eq x M]
        
        % Reshape Btdraw from [K x T] to [K_per_eq x M x T]
        Btdraw_reshaped = reshape(Btdraw, [K_per_eq, M, t]);
        
        % Permute to [T x K_per_eq x M] to match desired format
        Btdraw_4d = permute(Btdraw_reshaped, [3, 1, 2]);  % [T x K_per_eq x M]
        
        % Store results for this iteration
        residuals_save(counter,:,:) = Y_residuals';     % Transpose: [t x M]
        fitted_save(counter,:,:) = Y_fitted';           % Transpose: [t x M]
        Btdraw_save(counter,:,:,:) = Btdraw_4d;         % [t x K_per_eq x M]
        
    end % END saving after burn-in results 
end %END main Gibbs loop (for irep = 1:nrep+nburn)
clc;
toc; % Stop timer and print total time



end