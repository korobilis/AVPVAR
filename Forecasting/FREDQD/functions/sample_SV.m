function [sigma_t, isigma_t, h, sig] = sample_SV(SV_type, yhathat, h, sig, T, N)


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
        SSE       = diag(yhathat' * yhathat); % sum of squared errors per series
        r1        = (0.01 + T) / 2;
        r2        = (0.01 + SSE) / 2;
        iR        = gamrnd(r1, 1 ./ r2);       % inverse of variances
        isigma_t  = repmat(sqrt(iR'), T, 1);
        sigma_t   = 1 ./ isigma_t;
        h         = [];  % not used
        sig       = [];  % not used

    otherwise
        error('Unknown SV_type: %s. Use "SV" or "CV".', SV_type);
end

end
