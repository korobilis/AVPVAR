% TVP-VAR Time varying structural VAR with stochastic volatility
% ------------------------------------------------------------------------------------
% This code implements the TVP-VAR model as in Primiceri (2005). See also
% the monograph, Section 4.2 and Section 3.3.2 (Korobilis monograph)
% =========================================================================
% Written by: 
%
%    Dimitris Korobilis and modified by Nicolas Hardy
%    University of Glasgow
%
% First version: 20 April   2025
% This version:  20 April   2025
% =========================================================================
% ************************************************************************************
% The model is:
%
%     Y(t) = B0(t) + B1(t)xY(t-1) + B2(t)xY(t-2) + e(t) 
% 
%  with e(t) ~ N(0,SIGMA(t)), and  L(t)' x SIGMA(t) x L(t) = D(t)*D(t),
%             _                                          _
%            |    1         0        0       ...       0  |
%            |  L21(t)      1        0       ...       0  |
%    L(t) =  |  L31(t)     L32(t)    1       ...       0  |
%            |   ...        ...     ...      ...      ... |
%            |_ LN1(t)      ...     ...    LN(N-1)(t)  1 _|
% 
% 
% and D(t) = diag[exp(0.5 x h1(t)), .... ,exp(0.5 x hn(t))].
%
% The state equations are
%
%            B(t) = B(t-1) + u(t),            u(t) ~ N(0,Q)
%            l(t) = l(t-1) + zeta(t),      zeta(t) ~ N(0,S)
%            h(t) = h(t-1) + eta(t),        eta(t) ~ N(0,W)
%
% where B(t) = [B0(t),B1(t),B2(t)]', l(t)=[L21(t),...,LN(N-1)(t)]' and
% h(t) = [h1(t),...,hn(t)]'.
%
% ************************************************************************************
%   NOTE: 
%      There are references to equations of Primiceri, "Time Varying Structural Vector
%      Autoregressions & Monetary Policy",(2005),Review of Economic Studies 72,821-852
%      for your convenience. The definition of vectors/matrices is also based on this
%      paper.
% ------------------------------------------------------------------------------------

function [yfore_save] = TVP_RW_EB(Y,p,nrep,nburn,nfore)

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

%----------------------------PRELIMINARIES---------------------------------
% Set some Gibbs - related preliminaries
% nrep = 500;  % Number of replications
% nburn = 200;   % Number of burn-in-draws
% it_print = 100;  %Print in the screen every "it_print"-th iteration

%========= PRIORS:
% To set up training sample prior a-la Primiceri, use the following subroutine
[B_OLS,VB_OLS,A_OLS,sigma_OLS,VA_OLS]= ts_prior(Y,tau,M,p);

% % Or use uninformative values
% A_OLS = zeros(numa,1);
% B_OLS = zeros(K,1);
% VA_OLS = eye(numa);
% VB_OLS = eye(K);
% sigma_OLS = 0*ones(M,1);

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
    
    %% ----------------------------Forecasts AFTER-BURN-IN DRAWS  -----------------
    if irep > nburn 
        counter = counter + 1;
        [M, T] = size(y);
        B_T = Btdraw(:, T);
        H_T = Ht((T-1)*M+1:T*M, :);

        % Storage
        Y_fore = zeros(M, nfore);
        Ylags = y(:, T-p+1:T);

        for ifore = 1:nfore
              % Create regressor X_t = [1; y_{t-1}; ...; y_{t-p}]
              X_t = 1;
              for lag = 1:p
                 X_t = [X_t; Ylags(:, end - lag + 1)];
              end
          % Reshape B_T into matrix form: M x K
        Bmat = reshape(B_T, M, K/M);  % [M x K]

        % Forecast: y_t+1 = Bmat * X_t + e_t
        e_T = mvnrnd(zeros(M,1), H_T, 1)';  % M x 1
        y_next = Bmat * X_t + e_T;         % M x 1
        Y_fore(:, ifore) = y_next;

         % Update lags
         Ylags = [Ylags, y_next];
        end
         
       
      yfore_save(:,:,counter) = Y_fore;   
    end % END saving after burn-in results 
end %END main Gibbs loop (for irep = 1:nrep+nburn)
clc;
toc; % Stop timer and print total time






