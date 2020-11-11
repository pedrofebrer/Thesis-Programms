function V = Vega(K, S0, r , q, t, sig)

    % Inputs:
    % K - strike price vector
    % S0 - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % t - time to maturity value
    % sig - volatility value vector

    % Outputs:
    % V - Vega measure values vector

    d1 = (log(S0./K) + (r-q+(1/2)*sig.^2)*t)./(sig*sqrt(t));    % d1 is defined as a vector for the values of d_+ of the Black-Scholes formula for each value in K
    V = S0 .* exp(-q.*t) .* normpdf(d1).* sqrt(t) ;                 % V is defined as a vector for the values of the Vega nmeasure under the Black-Scholes model for each value in K 
end