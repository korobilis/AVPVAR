function [beta_save,sigma_save,Gamma_save] = AP_RW(y,X,Z,ngibbs,nburn)

[T,p] = size(X);
%Z = normalize(Z);
M = size(Z,2);
Mp = (M+1)*p;

algo = 1 + double(Mp<T);

% Some useful matrices for the time-varying parameters
Tp = T*p;
x = []; xZ = [];
H         = speye(T,T) - sparse(2:T,1:(T-1),ones(1,(T-1)),T,T);
Hinv      = speye(T,T)/H;
for i = 1:p
    SURx = full(SURform(X(:,i))*Hinv);
    x  = [x, SURx];
    xZ = [xZ, SURx*Z];
end


% Priors for Gamma
Q    = .01*ones(Mp,1);

% Horseshoe prior
lambda = ones(Mp,1);
tau = 1;
nu = ones(Mp,1);
xi = 1;
In = eye(T);

% initialize volatility parameters
isigma_t = 0.1*ones(T,1);
h        = ones(T,1); 
sig      = 0.1;

%% ========|STORAGE MATRICES:
beta_save  = zeros(ngibbs,T,p);
Gamma_save = zeros(ngibbs,Mp);
sigma_save = zeros(ngibbs,T);
%% ======================= BEGIN MCMC ESTIMATION =======================
disp('Run Gibbs sampler for InstrumentedTVP algorithm');
for iter = 1: (ngibbs + nburn)
    
    %% Sample gamma
    xs = isigma_t.*[X, xZ];
    ys = isigma_t.*y;
    Gamma = randn_gibbs(ys,xs,Q,T,Mp,In,algo);
    beta = zeros(T,p);
for i = 1:p
    beta(:,i) = Gamma(i) + cumsum(Z*Gamma((i-1)*M+p+1:i*M+p));
end


    %% sample Q
    [Q,~,tau,xi,nu] = horseshoe_prior(Gamma,Mp,tau,xi,nu);    

    %% Sample sigma2_t
    yhat = y - [X xZ]*Gamma;
    ystar  = log(yhat.^2 + 1e-6);        
    [h, ~] = SVRW(ystar,h,sig,4);  % log stochastic volatility
    isigma_t = exp(-0.5*h);
    sigma_t  = exp(h);
    r1 = 1 + T - 1;   r2 = 0.01 + sum(diff(h).^2)';
    sig = 1./gamrnd(r1./2,2./r2);   % sample state variance of log(sigma_t.^2)
         
    
    %% Save draws
    if iter > nburn
        beta_save(iter-nburn,:,:) = beta;
        Gamma_save(iter-nburn,:,:) = Gamma;
        sigma_save(iter-nburn,:) = sigma_t;
    end   
end
end
