    % Btdrawc is a draw of the mean VAR coefficients, B(t)
    [Btdrawc,log_lik] = carter_kohn(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar);
    
    % Accept draw
    Btdraw = Btdrawc;


    %=====| Draw Q, the covariance of B(t) (from iWishart)
    % Take the SSE in the state equation of B(t)
    Btemp = Btdraw(:,2:t)' - Btdraw(:,1:t-1)';
    sse_2 = zeros(K,K);
    for i = 1:t-1
        sse_2 = sse_2 + Btemp(i,:)'*Btemp(i,:);
    end
    
    % ...and subsequently draw Q, the covariance matrix of B(t)
    Qinv = inv(sse_2 + Q_prmean);
    Qinvdraw = wish(Qinv,t-1+Q_prvar);
    Qdraw = inv(Qinvdraw);  % this is a draw from Q