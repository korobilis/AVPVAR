function [htilde, h0, omegah, omegahhat, Domegah] = SVRW_gam_Nico(Ystar, htilde, h0, omegah, b0, Vh0, Vh, Vomegah)
% SVRW_GAM_NICO - Stochastic Volatility Random Walk Sampler
% ================================================================================
% MODEL DESCRIPTION
% ================================================================================
% This function implements the auxiliary mixture sampler for stochastic volatility
% models, following Kim, Shephard, and Chib (1998) and Omori et al. (2007).
%
% 1. ORIGINAL SV MODEL:
%    y*_t = exp(h_t/2) * ε_t,  where ε_t ~ N(0,1)
%    
%    Taking logs and squaring:
%    log(y*_t^2) = h_t + log(ε_t^2)
%    
%    where log(ε_t^2) ~ log(χ²_1) (log chi-squared with 1 df)
%
% 2. MIXTURE APPROXIMATION:
%    The distribution log(χ²_1) is approximated by a 7-component normal mixture:
%    log(ε_t^2) ≈ Σ_{j=1}^7 p_j * N(m_j - 1.2704, σ_j^2)
%    
%    where p_j are mixture probabilities, m_j are means, σ_j^2 are variances
%    The constant 1.2704 = E[log(χ²_1)] ensures correct mean
%
% 3. HIERARCHICAL REPRESENTATION:
%    Y*_t = h_0 + ω_h * h̃_t + m_{S_t} + ν_t
%    
%    where:
%    - Y*_t = log(y*_t^2 + c), with small constant c for numerical stability
%    - h_t = h_0 + ω_h * h̃_t (location-scale representation)
%    - S_t ∈ {1,...,7} is the mixture component indicator
%    - ν_t ~ N(0, σ_{S_t}^2)
%    
% 4. STATE EQUATION FOR h̃_t:
%    h̃_t = h̃_{t-1} + u_t,  where u_t ~ N(0,1)
%    
%    In matrix form: H*h̃ = ν, where ν ~ N(0, S)
%    - H is the difference matrix
%    - S = diag(V_h, 1, ..., 1)
%
% 5. PRIORS:
%    - h_0 ~ N(b_0, V_{h0})
%    - ω_h ~ N(0, V_{ωh})
%    - h̃_1 ~ N(0, V_h)
%
% ================================================================================
% INPUTS:
%   Ystar   - T×1 vector of transformed observations log(y*^2 + c)
%   htilde  - T×1 current draw of standardized volatility process
%   h0      - Scalar, current draw of initial volatility level
%   omegah  - Scalar, current draw of volatility scaling parameter
%   b0      - Prior mean for h_0
%   Vh0     - Prior variance for h_0
%   Vh      - Prior variance for h̃_1
%   Vomegah - Prior variance for ω_h
%
% OUTPUTS:
%   htilde   - T×1 new draw of standardized volatility process
%   h0       - New draw of initial volatility level
%   omegah   - New draw of volatility scaling parameter
%   omegahhat - Posterior mean of ω_h (for density estimation)
%   Domegah  - Posterior variance of ω_h
% ================================================================================

    %% ============================================================
    %% STEP 1: SETUP AND MIXTURE PARAMETERS
    %% ============================================================
    
    T = length(htilde);  % Number of time periods
    
    % Kim-Shephard-Chib 7-component mixture parameters
    % These approximate the log(χ²_1) distribution
    
    % Mixture probabilities (sum to 1)
    pj = [0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.2575];
    % Mixture means (already adjusted by -1.2704 to have correct expectation)
    mj = [-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819] - 1.2704;
    % Mixture variances
    sigj = [5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261];
    % Pre-compute standard deviations for efficiency
    sqrtsigj = sqrt(sigj);
    
    %% ============================================================
    %% STEP 2: SAMPLE MIXTURE INDICATORS S_t
    %% ============================================================
    % For each t, sample S_t from discrete distribution:
    % P(S_t = j | Y*_t, h_t) ∝ p_j * N(Y*_t | h_t + m_j, σ_j^2) 
    % Generate uniform random numbers for inverse CDF method
    temprand = rand(T, 1);
    % Compute unnormalized probabilities for each component j and time t
    q = repmat(pj, T, 1) .* ...
        normpdf(repmat(Ystar, 1, 7), ...                     % Data
                repmat(h0 + omegah*htilde, 1, 7) + ...       % Mean: h_t + m_j
                repmat(mj, T, 1), ...
                repmat(sqrtsigj, T, 1));                     % Std dev: σ_j
    
    % Normalize to get probabilities (each row sums to 1)
    q = q ./ repmat(sum(q, 2), 1, 7);
    
    % Sample S_t using inverse CDF method
    S = 7 - sum(repmat(temprand, 1, 7) < cumsum(q, 2), 2) + 1;
    
    %% ============================================================
    %% STEP 3: SAMPLE h̃ CONDITIONAL ON S, h_0, ω_h
    %% ============================================================
    % Model after conditioning on S:
    % Y*_t = h_0 + ω_h*h̃_t + d_t + ε_t
    
    % Setup precision matrices
    Hh = speye(T) - sparse(2:T, 1:(T-1), ones(1,T-1), T, T);
    invSh = sparse(1:T, 1:T, [1/Vh; ones(T-1,1)]);
    
    % Setup observation equation parameters
    dconst = mj(S)';  % d_t = m_{S_t}
    invOmega = sparse(1:T, 1:T, 1./sigj(S));  % Diagonal precision matrix
    
    % Compute posterior precision and mean
    Kh = Hh'*invSh*Hh + invOmega*omegah^2;
    htildehat = Kh\(invOmega*omegah*(Ystar - dconst - h0));
    
    % Sample from posterior
    htilde = htildehat + chol(Kh,'lower')'\randn(T,1);
    
    %% ============================================================
    %% STEP 4: SAMPLE (h_0, ω_h) JOINTLY
    %% ============================================================
    % Linear regression: Y*_t - d_t = h_0 + ω_h*h̃_t + ε_t
    
    % Setup design matrix and prior
    Xbeta = [ones(T,1), htilde];
    invVbeta = diag([1/Vh0, 1/Vomegah]);
    
    % Compute posterior parameters
    XbetainvOmega = Xbeta'*invOmega;
    invDbeta = invVbeta + XbetainvOmega*Xbeta;
    betahat = invDbeta\(invVbeta*[b0; 0] + XbetainvOmega*(Ystar - dconst));
    
    % Sample from posterior
    beta = betahat + chol(invDbeta,'lower')'\randn(2,1);
    h0 = beta(1);
    omegah = beta(2);
    
    %% ============================================================
    %% STEP 5: SIGN NORMALIZATION
    %% ============================================================
    % The model h_t = h_0 + ω_h*h̃_t is invariant to sign changes
    % We break this symmetry by randomly flipping signs with prob 0.5
    
    U = -1 + 2*(rand > 0.5);  % Random sign: +1 or -1
    htilde = U*htilde;         % Flip h̃
    omegah = U*omegah;        % Flip ω_h
    
    %% ============================================================
    %% STEP 6: COMPUTE POSTERIOR MOMENTS FOR ω_h
    %% ============================================================
    % These are useful for posterior density estimation
    
    % Re-compute posterior covariance matrix with sign-corrected values
    Xbeta = [ones(T,1), htilde];
    Dbeta = (invVbeta + Xbeta'*invOmega*Xbeta)\speye(2);
    
    % Extract marginal posterior moments for ω_h
    omegahhat = betahat(2);    % Posterior mean of ω_h
    Domegah = Dbeta(2,2);      % Posterior variance of ω_h
    
end