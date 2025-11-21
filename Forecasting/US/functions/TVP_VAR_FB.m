function  [yfore_save,beta_save,invalid_runs,beta_OLS,residuals_save] = TVP_VAR_FB(Y,nfore,p,r,ngibbs,nburn,nthin)
%% =========================================================================
% Written by: 
%
%    Dimitris Korobilis and Nicolas Hardy
%    University of Glasgow
%
% First version: 05 May 2025
% This version : 08 Nov 2025
%% =========================================================================
%=========| Estimation
[y,X,N,T,~,K] = prepare_BVAR_matrices(Y,p,0);
        
%% ==================| 1. OLS estimates:
[beta_OLS,~,~] = getOLS(y,X,N,p);
%==================| 2. PRIORS:
In = eye(T);
%% Extract factors
yhat       = y - X*beta_OLS;
[F,L]      = extract(zscore(yhat),r);
F          = F/chol(cov(F)) - mean(F/chol(cov(F)));
L          = repmat(L,1,1,T);
L          = permute(L,[3 2 1]); % L(T,r,N)
%% Volatility parameters (Initialise)
isigma_t   = 0.1*ones(T,N);
sigma_t    = 1./isigma_t;
h          = ones(T,N); 
sig        = 0.1.*ones(1,N);
%% Hyhperparameters for Loadings-Instrumented equation ~ Horseshoe
psi_L = cell(N, 1);
for ieq = 1:N
    psi_L{ieq} = ones(r * T, 1);
end
tau_L = ones(N, 1);
%% ========== Prepare matrices for TVPs and TVLs
H_L    = speye(T*r) - sparse(r+1:T*r,1:(T-1)*r,1,(T*r),(T*r));
Hinv_L = speye(T*r)/H_L;
x  = []; H = []; Hinv = [];
TK = T*K;
H    = speye(TK) - sparse(K+1:TK,1:(T-1)*K,1,TK,TK);
Hinv = speye(TK)/H;
x    = full(SURform(X) * Hinv);

%% ========== Some priors 
MK   = T * K;
psi  = cell(N, 1);
for ieq = 1:N
    psi{ieq} = ones(MK, 1);
end
tau  = ones(N, 1);
beta = zeros(T, K, N);
%========|STORAGE MATRICES:
beta_save  = zeros(ngibbs/nthin,T,K,N);
L_save     = zeros(ngibbs/nthin,T,r,N); 
F_save     = zeros(ngibbs/nthin,T,r);
R_save     = zeros(ngibbs/nthin,T,N);
invalid_runs = 0;

%% ===================== Begin Gibbs Sampler ==========================
tic;
counter = 0;
fprintf('\n')
fprintf('TVP-VAR-FB horseshoe');
fprintf('Iteration 000000')
for iter = 1:ngibbs + nburn
    % Print every 100th iterations on the screen
    if mod(iter,100) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%6d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',iter)        
    end

    %% ==================== Sample coefficients Beta =======================
    for ieq = 1:N
        for t = 1:T
            y_til(t,ieq) = y(t,ieq) - L(t,:,ieq)*F(t,:)';  %  L size (T x r x N)
        end
    end

    stationary = 0;
    count = 0;
    invalid = 0;

    while ~stationary
        beta = zeros(T, K, N); % reinitialize beta
    
        for ieq = 1:N
        % STEP 1.1: Construct design matrix xs
            xs = isigma_t(:, ieq) .* x;
        % STEP 1.2: Construct dependent variable
            ys = isigma_t(:, ieq) .* y_til(:, ieq);
        % STEP 1.3: Sample Gamma using Horseshoe prior
            [Gamma(:,ieq), psi{ieq}, tau(ieq)] = horseshoeG(ys, xs, In, psi{ieq}, tau(ieq),1, 1 + (size(xs,1) > size(xs,2)));                                    
        % STEP 1.4: Recover beta
            beta(:, :, ieq) = reshape(Hinv * Gamma(:,ieq), K, T)';
        end
    
        % STEP 2: Stationarity check — applied to ALL TVP types
        B = [squeeze(beta(T, 2:end, :))'; eye(K - 1 - N), zeros(K - 1 - N, N)];
    
        if max(abs(eig(B))) < 0.999
            stationary = 1;
        else
            count = count + 1;
            if count > 4000
                % Multiply Gamma by 0.999 to push toward stationarity
                fprintf('Warning: Stationarity not achieved after 4000 tries. Scaling Gamma by 0.999 and continuing...\n');
                Gamma = Gamma * 0.999;
                invalid = 1;
                % Re-derive beta from the modified Gamma
                for ieq = 1:N
                    beta(:, :, ieq) = reshape(Hinv * Gamma(:,ieq), K, T)';
                end
                % Reset count and continue the loop to check stationarity again
                count = 0;
            end
        end
    end
    invalid_runs = invalid_runs + invalid;
      
    %% ======== Step 2: Sample Factors and Loadings =========
    for ieq = 1:N  
        yhat(:,ieq) = y(:,ieq) - dot(X,beta(:,:,ieq),2);  
    end

    %% ============ Sample Factors   ==========
    for t = 1:T
        LL = reshape(L(t,:,:), r, N);  
        F_var  = eye(r)/( eye(r) + LL*diag(isigma_t(t,:).^2)*LL' );
        F_mean = F_var * LL*diag(isigma_t(t,:).^2)*yhat(t,:)';
        F(t,:) = (F_mean + chol(F_var)'*randn(r,1))';
    end
    F = (F-mean(F))/chol(cov(F));

    %% ============ Sample Loadings   ==========
    f  = []; 
    f = full(SURform(F) * Hinv_L);
    % ============ Step 2: Sample Delta and recover L ============
    L = zeros(T, r, N);
    for ieq = 1:N
    % Step 2.1: Construct fZs (design matrix)
        fs = isigma_t(:,ieq) .* f;
    % Step 2.2: Construct dependent variable
        ys2 = isigma_t(:,ieq) .* yhat(:,ieq);
    % Step 2.3: Sample Delta using Horseshoe prior
        [Delta(:,ieq), psi_L{ieq}, tau_L(ieq)] = horseshoeG(ys2, fs, In,psi_L{ieq}, tau_L(ieq), 1, 1 + (size(fs,1) > size(fs,2)));
    % Step 2.4: Recover Loadings
        L(:,:,ieq) = reshape(Hinv_L * Delta(:,ieq), r, T)';
    end
    %% ========== Stochastic Volatility Block======
    yhathat = 0*y;
    for ieq = 1:N
        for t = 1:T
            yhathat(t,ieq) = yhat(t,ieq) - L(t,:,ieq)*F(t,:)';  % L(TxrxN)
        end   
    end      
    ystar = log(yhathat.^2 + 1e-10);
    for ieq = 1:N
        h(:, ieq) = SVRW(ystar(:, ieq), h(:, ieq), sig(:, ieq), 1);
        sigma_t(:, ieq)  = exp(0.5 * h(:, ieq));
        isigma_t(:, ieq) = exp(-0.5 * h(:, ieq));
    % Sample state variance of log-volatility innovations
        r1 = 1 + T - 1;
        r2 = 0.01 + sum(diff(h(:, ieq)).^2);
        sig(:, ieq) = 1 ./ gamrnd(r1 / 2, 2 / r2);
    end
       
        %% ============ Save and Forecasting ==================  
        if iter>nburn && mod(iter,nthin) == 0 
            counter = counter + 1;
            beta_save(counter,:,:,:)    = beta; 
            L_save(counter,:,:,:)       = L;    
            F_save(counter,:,:)         = F;
            R_save(counter,:,:)         = sigma_t;
            residuals_save(counter,:,:) = yhathat;

            % ============| FORECASTING |============       
            beta_mat = squeeze(beta(end,:,:));
            LL = reshape(L(end,:,:), r, N);  
            SIGMAf = LL'*(F(end,:)'*F(end,:))*LL + diag(sigma_t(end,:).^2);  
            cSIGMA = chol(SIGMAf);

            y_fore = zeros(nfore,N);        
            x_fore = [y(end,:) X(end,2:end-N)];  
            y_fore(1,:) = [1 x_fore]*beta_mat + randn(1,N)*cSIGMA;     
            for ifore = 2:nfore           
                x_fore          = [y_fore(ifore-1,:) x_fore(:,1:end-N)];        
                y_fore(ifore,:) = [1 x_fore]*beta_mat  + randn(1,N)*cSIGMA;   
            end        
            yfore_save(:,:,counter) = y_fore';
        end
end
toc
end