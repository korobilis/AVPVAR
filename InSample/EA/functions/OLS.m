function [sigma,se,beta] = OLS(X,Y)

beta = (X'*X)\(X'*Y);
sigma = (Y-X*beta)'*(Y-X*beta)/(size(Y,1)-size(X,2));
se = sqrt(diag(sigma*inv(X'*X)));