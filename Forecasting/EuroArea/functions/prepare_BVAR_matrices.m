function [Y,X,N,T,KK,K] = prepare_BVAR_matrices(Y,p,tvp)

%% Outputs
% T: number of observations in the times series Y
% N: number of variables in the VAR

% Number of observations and dimension of X and Y
T = size(Y,1); % T is the time-series observations of Y
N = size(Y,2); % N is the dimensionality of Y

if tvp
    K = N*p + T-p;
else
    K = N*p + 1;
end
KK = K*N;
	
% ===================================| VAR EQUATION |==============================
% Take lags, and correct observations
Ylag = mlag2(Y,p);
Ylag = Ylag(p+1:end,:);
T = T-p;

% Final VAR matrices/vectors
% 1/ VAR matrices for traditional matrix form
if tvp
    X = [tril(ones(T)) Ylag];
else
    X = [ones(T,1) Ylag];
end
Y = Y(p+1:end,:);