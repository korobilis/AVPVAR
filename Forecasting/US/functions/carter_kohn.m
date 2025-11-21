function [bdraw,log_lik] = carter_kohn(y,Z,Ht,Qt,m,p,t,B0,V0)
%% Carter and Kohn (1994), On Gibbs sampling for state space models.
% Y_t = Z_t*theta_t + e_t; V(et) = H_t
% Standard KF notation: Y_t = H*theta_t + e_t; V(et) = R
%% Kalman Filter
bp = B0;
Vp = V0;
bt = zeros(t,m);
Vt = zeros(m^2,t);
log_lik = 0;
for i=1:t
    R = Ht((i-1)*p+1:i*p,:);  % variance measurement equation
    H = Z((i-1)*p+1:i*p,:);
    cfe = y(:,i) - H*bp;      % conditional forecast error
    f = H*Vp*H' + R;          % variance of the conditional forecast error
    K = Vp*H'*inv(f);          % Kalman Gain
    btt = bp + K*cfe;         % Update b_(t|t)
    Vtt = Vp-K*H*Vp;         %Update P_(t|t)
    if i < t
        bp = btt;             % Prediction step 
        Vp = Vtt + Qt;        % State prediction error variance
    end
    bt(i,:) = btt';
    Vt(:,i) = reshape(Vtt,m^2,1);
end

%% draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw = zeros(t,m);
bdraw(t,:) = mvnrnd(btt,Vtt,1);

%% Backward recurssions
for i=1:t-1
    bf = bdraw(t-i+1,:)';
    btt = bt(t-i,:)';
    Vtt = reshape(Vt(:,t-i),m,m);
    f = Vtt + Qt;
    inv_f = inv(f);
    cfe = bf - btt;
    bmean = btt + Vtt*inv_f*cfe;
    bvar = Vtt - Vtt*inv_f*Vtt;
    bdraw(t-i,:) = mvnrnd(bmean,bvar,1); 
end
bdraw = bdraw';