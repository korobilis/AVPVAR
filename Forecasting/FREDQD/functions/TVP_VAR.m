function  [yfore_save,beta_save,F_save,L_save,R_save,beta_OLS] = TVP_VAR(Y,nfore,p,r,ngibbs,nburn,nthin,tvp_type,tvp_loadings,SV_type)
%% =========================================================================
% This function allows to estimate a TVP model with stochastic volatility
% The configurations are: 
% tvp_type = 'TVP-RW' or 'CP' (Parameters evolving as random walks, or Constant Parameters)
% tvp_type = 'TVL-RW' or 'CL' (Loadings evolving as random walks, or Constant Loadings)
% SV_type = 'SV' or 'CV' (Stochastic Volatility, or Constant volatility
% Written by: 
%
%    Dimitris Korobilis and Nicolas Hardy
%    University of Glasgow
%
% First version: 05 May 2025
% This version : 18 May 2025
%% =========================================================================
%=========| Estimation
[y,X,N,T,~,K] = prepare_BVAR_matrices(Y,p,0);
%% ==================| 1. OLS estimates:
[beta_OLS,~,~] = getOLS(y,X,N,p);
%% ==================| 2. PRIORS:
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
switch tvp_loadings
    case 'CL'
        for ieq = 1:N
            psi_L{ieq} = ones(r, 1);
        end
        tau_L = ones(N, 1);
    case 'TVL-RW'
        for ieq = 1:N
            psi_L{ieq} = ones(r * T, 1);
        end
        tau_L = ones(N, 1);
end
%% ========== Prepare matrices for TVPs and TVLs
switch tvp_loadings
    case 'TVL-RW'
        H_L    = speye(T*r) - sparse(r+1:T*r,1:(T-1)*r,1,(T*r),(T*r));
        Hinv_L = speye(T*r)/H_L;

    case 'CL'
        H_L    = [];
        Hinv_L = [];
end
% ----------- TVP type case -----------
x  = [];xZ = [];H = [];Hinv = [];TK = T*K;
switch tvp_type
    case 'CP'
        for i = 1:K
            SURx_i = full(SURform(X(:,i)));
            x      = [x, SURx_i];
            xZ     = [xZ, SURx_i * ones(T,1)];
        end
    case 'TVP-RW'
        H    = speye(TK) - sparse(K+1:TK,1:(T-1)*K,1,TK,TK);
        Hinv = speye(TK)/H;
        x    = full(SURform(X) * Hinv);
end
%% ========== Some priors: Autorregressive Parameter Block
switch tvp_type
    case 'TVP-RW'
        MK   = T * K;
        psi  = cell(N, 1);
        for ieq = 1:N
            psi{ieq} = ones(MK, 1);
        end
        tau  = ones(N, 1);
        beta = zeros(T, K, N);
    case 'CP'
        psi  = cell(N, 1);
        for ieq = 1:N
            psi{ieq} = ones(K, 1);
        end
        tau  = ones(N, 1);
        beta = zeros(T, K, N);
end

%========|STORAGE MATRICES:
beta_save  = zeros(ngibbs/nthin,T,K,N);
L_save     = zeros(ngibbs/nthin,T,r,N); 
F_save     = zeros(ngibbs/nthin,T,r);
R_save     = zeros(ngibbs/nthin,T,N);

%% ===================== Begin Gibbs Sampler ==========================
tic;
counter = 0;
fprintf('\n')
fprintf('TVP-VAR with %s, %s and %s\n', char(tvp_type), char(tvp_loadings), char(SV_type));
fprintf('Iteration 000000')
for iter = 1:ngibbs + nburn
    % Print every 100th iterations on the screen
    if mod(iter,100) == 0
        fprintf('%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%s%6d',8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,'Iteration ',iter)
        
    end

    %% ==================== Sample coefficients Beta =======================
        for ieq = 1:N
            for t=1:T
                y_til(t,ieq) = y(t,ieq) - L(t,:,ieq)*F(t,:)';  %  L size (T x r x N)
            end
        end

    stationary = 0;
    count = 0;
    invalid = 0;

    while ~stationary
        beta = zeros(T, K, N);  % reinitialize beta
        for ieq = 1:N
        % STEP 1.1: Construct design matrix xs
            switch tvp_type
                case 'TVP-RW'
                    xs = isigma_t(:, ieq) .* x;
                case 'CP'
                    xs = isigma_t(:, ieq) .* xZ;
            end
        % STEP 1.2: Construct dependent variable
            ys = isigma_t(:, ieq) .* y_til(:, ieq);
        % STEP 1.3: Sample Gamma using Horseshoe prior
            [Gamma(:,ieq), psi{ieq}, tau(ieq)] = horseshoeG(ys, xs, In, psi{ieq}, tau(ieq), ...
                1, 1 + (size(xs,1) > size(xs,2)));
        % STEP 1.4: Recover beta
            switch tvp_type
                case 'CP'
                    for i = 1:K
                        beta(:, i, ieq) = Gamma(i,ieq);
                    end
                case 'TVP-RW'
                    beta(:, :, ieq) = reshape(Hinv * Gamma(:,ieq), K, T)';
            end
        end

    % STEP 2: Stationarity check — applied to ALL TVP types
        B = [squeeze(beta(T, 2:end, :))'; eye(K - 1 - N), zeros(K - 1 - N, N)];
        if max(abs(eig(B))) < 0.999
            stationary = 1;
        else
            count = count + 1;
            if count > 4000
                warning('Stationarity not achieved after 4000 tries. Continuing with last draw.');
                B = B.*0.99;
                invalid=invalid+1;
                break;
            end
        end
    end
      
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
        f  = []; fZ = [];

        switch tvp_loadings
            case 'CL'
                for i = 1:r
                    SURf = full(SURform(F(:,i)));
                    f  = [f, SURf];
                    fZ = [fZ, SURf * ones(T,1)];
                end
            case 'TVL-RW'
                f = full(SURform(F) * Hinv_L);
        end

        % ============ Sample Delta and recover L ============
        L = zeros(T, r, N);

        for ieq = 1:N
    % Step 2.1: Construct fZs (design matrix)
            switch tvp_loadings
                case 'CL'
                    fZs = isigma_t(:,ieq) .* fZ;
                case 'TVL-RW'
                    fZs = isigma_t(:,ieq) .* f;
            end

    % Step 2.2: Construct dependent variable
            ys2 = isigma_t(:,ieq) .* yhat(:,ieq);
    % Step 2.3: Sample Delta using Horseshoe prior
            [Delta(:,ieq), psi_L{ieq}, tau_L(ieq)] = horseshoeG(ys2, fZs, In, ...
                psi_L{ieq}, tau_L(ieq), 1, 1 + (size(fZs,1) > size(fZs,2)));
    % Step 2.4: Recover Loadings
            switch tvp_loadings
                case 'CL'
                    for i = 1:r
                        L(:,i,ieq) = Delta(i,ieq);
                    end
                case 'TVL-RW'
                    L(:,:,ieq) = reshape(Hinv_L * Delta(:,ieq), r, T)';
            end
        end

        %% ========== Stochastic Volatility Block======
        yhathat = 0*y;
        for ieq = 1:N
            for t=1:T
                yhathat(t,ieq) = yhat(t,ieq) - L(t,:,ieq)*F(t,:)';  
            end   
        end   

        switch SV_type
            case 'SV'
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

            case 'CV'
                SSE       = diag(yhathat' * yhathat); 
                r1        = (0.01 + T) / 2;
                r2        = (0.01 + SSE) / 2;
                iR        = gamrnd(r1, 1 ./ r2);       
                isigma_t  = repmat(sqrt(iR'), T, 1);
                sigma_t   = 1 ./ isigma_t;
        end
       
        %% ============ Save and Forecasting ==================  
        if iter>nburn && mod(iter,nthin) == 0 
            counter = counter + 1;
            beta_save(counter,:,:,:)  = beta; 
            L_save(counter,:,:,:)     = L;    
            F_save(counter,:,:)       = F;
            R_save(counter,:,:)       = sigma_t;
            % ------------------| FORECASTING |------------------        
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