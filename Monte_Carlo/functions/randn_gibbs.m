function [Beta] = randn_gibbs(y,X,lambda,n,p,I_n,algorithm)

if algorithm == 1
    U = bsxfun(@times,lambda,X');
    %% step 1 %%
    u = normrnd(zeros(p,1),sqrt(lambda));	
    v = X*u + randn(n,1);
    %% step 2 %%
    v_star = ((X*U) + I_n)\(y-v);
    Beta = (u + U*v_star);
elseif algorithm == 2
    %% matrices %%
    Q_star = X'*X;
    Dinv = diag(1./lambda);       
    L=chol((Q_star + Dinv),'lower');
    v=L\(y'*X)';
    mu=L'\v;
    u=L'\randn(p,1);
    Beta = mu+u;
end