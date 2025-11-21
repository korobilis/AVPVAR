function  [yfore_save,beta_save,invalid_runs,beta_OLS,residuals_save,Gamma_save,tau_save,psi_save] = APVAR(Y,Z,nfore,p,r,ngibbs,nburn,nthin,SV_type)
%% =========================================================================
% Written by: 
%
%    Dimitris Korobilis and Nicolas Hardy
%    University of Glasgow
%
% First version: 05 May 2025
% This version : 07 November 2025
%% =========================================================================

%% ===================| Prepare data |=============================================
[y,X,N,T,~,K] = prepare_BVAR_matrices(Y,p,0);
Zfore = Z(p+1:end,:);
[Z,M]         = prepare_Z2(Z,p) ; 
ML = size(Z,2);

% ====================| 1. OLS estimates:
[beta_OLS,~,~] = getOLS(y,X,N,p);
% ====================| 2. PRIORS and EXTRACT FACTORS:
In = eye(T);
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
%% Hyhperparameters for Adaptive Loadings equation ~ Horseshoe
psi_L = cell(N, 1);
     for ieq = 1:N
         psi_L{ieq} = ones(r * (ML + 1), 1);
     end
tau_L = ones(N, 1);

%% ========== Prepare matrices for TVPs and TVLs
H_L    = speye(T) - sparse(2:T,1:T-1,1,T,T);
Hinv_L = speye(T)/H_L;
x  = [];xZ = [];H = [];Hinv = [];
TK = T*K;
H    = speye(T) - sparse(2:T,1:T-1,1,T,T);
Hinv = speye(T)/H;
     for i = 1:K
         SURx_i = full(SURform(X(:,i)) * Hinv);
         x      = [x, SURx_i];
         xZ     = [xZ, SURx_i * Z];
     end
%% Hyhperparameters for Adaptive Parameters equation ~ Horseshoe
MK   = (M + 1) * K;
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
fprintf('AP-VAR horseshoe with %s\n', char(SV_type));
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
            beta = zeros(T, K, N); 
    
            for ieq = 1:N
                % STEP 1.1: Construct design matrix xs
                xs = isigma_t(:, ieq) .* [X, xZ];    
                % STEP 1.2: Construct dependent variable
                ys = isigma_t(:, ieq) .* y_til(:, ieq);       
                % STEP 1.3: Sample Gamma using Horseshoe prior
                [Gamma(:,ieq), psi{ieq}, tau(ieq)] = horseshoeG(ys, xs, In, psi{ieq}, tau(ieq), ...
                                                          1, 1 + (size(xs,1) > size(xs,2)));       
                % STEP 1.4: Recover beta
                for i = 1:K
                    beta(:, i, ieq) = Gamma(i,ieq) + cumsum(Z * Gamma((i-1)*M + K + 1 : i*M + K,ieq));
                end
            end  

            % STEP 2: Stationarity check 
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
                                for i = 1:K
                                    beta(:, i, ieq) = Gamma(i,ieq) + cumsum(Z * Gamma((i-1)*M + K + 1 : i*M + K,ieq));
                                end
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

        %% ============ Step 2.1: Sample Factors   ==========
        for t = 1:T
            LL = reshape(L(t,:,:), r, N);  
            F_var  = eye(r)/( eye(r) + LL*diag(isigma_t(t,:).^2)*LL' );
            F_mean = F_var * LL*diag(isigma_t(t,:).^2)*yhat(t,:)';
            F(t,:) = (F_mean + chol(F_var)'*randn(r,1))';
        end
        F = (F-mean(F))/chol(cov(F));

        %% ============ | Step 2.1: Sample Loadings | =============
        f  = []; fZ = [];

        for i = 1:r
            SURf = full(SURform(F(:,i)) * Hinv_L);
            f  = [f, SURf];
            fZ = [fZ, SURf * Z];
        end
        % ============ | Step 2.2: Sample Delta and recover L | =============
        L = zeros(T, r, N);

        for ieq = 1:N
            %  Construct fZs (design matrix)
            fZs = isigma_t(:,ieq) .* [F, fZ];
            %  Construct dependent variable
            ys2 = isigma_t(:,ieq) .* yhat(:,ieq);
            % Sample Delta using Horseshoe prior
            [Delta(:,ieq), psi_L{ieq}, tau_L(ieq)] = horseshoeG(ys2, fZs, In, ...
                 psi_L{ieq}, tau_L(ieq), 1, 1 + (size(fZs,1) > size(fZs,2)));
            %  Recover Loadings
            for i = 1:r
                  L(:,i,ieq) = Delta(i,ieq) + cumsum(Z * Delta((i-1)*ML + r + 1 : i*ML + r,ieq));
            end
        end

        %% ========== Stochastic/Constant Volatility Block======
        yhathat = 0*y;
        for ieq = 1:N
            for t=1:T
                yhathat(t,ieq) = yhat(t,ieq) - L(t,:,ieq)*F(t,:)'; 
            end   
        end      
        [sigma_t, isigma_t, h, sig] = sample_SV(SV_type, yhathat, h, sig, T, N);
       
        %% ============ Save and Forecasting ==================  
        if iter>nburn && mod(iter,nthin) == 0 
            counter = counter + 1;
            beta_save(counter,:,:,:)  = beta; 
            L_save(counter,:,:,:)     = L;   
            F_save(counter,:,:)       = F;
            R_save(counter,:,:)       = sigma_t;
            residuals_save(counter,:,:) =yhathat;
            Gamma_save(counter,:,:) = Gamma; 
            Delta_save(counter,:,:) = Delta; 
            tau_save(counter,:)       = tau;
            tau_save_L(counter,:)     = tau_L;
            psi_save(counter,:) = psi;
         %% ============| FORECASTING |============     
           for ieq = 1:N
                for i = 1:K
                    temp = cumsum(Zfore * Gamma((i-1)*M + K + 1 : i*M + K,ieq));
                    beta_mat(i, ieq) = Gamma(i,ieq) + temp(end,1);
                end
           end

           for ieq = 1:N
                for i = 1:r
                    temp = cumsum(Zfore * Delta((i-1)*ML + r + 1 : i*ML + r,ieq));
                    Lend(i,ieq) = Delta(i,ieq) + temp(end,1);
                end
           end


              F_var_fore  = eye(r)/( eye(r) + Lend*diag(isigma_t(end,:).^2)*Lend' );
              F_mean_fore = F_var_fore * Lend*diag(isigma_t(end,:).^2)*yhat(end,:)';
              F_fore(1,:) = (F_mean_fore + chol(F_var_fore)'*randn(r,1))';
        
              LL = reshape(Lend, r, N);  
              SIGMAf = LL'*F_fore'*F_fore*LL + diag(sigma_t(end,:).^2);  
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