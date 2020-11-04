function sig = ImpliedVolatility (C, K, S0, r , q, t, sig0, tol, maxIterations)
    
    % Inputs:
    % C - European Call Option value vector
    % K - strike price vector
    % S - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % t - time to maturity value
    % sig0 - initial volatility value for the Newton Method
    % tol - tolerance error
    % maxIterations - maximum number of iterations of the Newton Method

    % Outputs:
    % sig - volatility values vector

    n = length(C);                  % n is defined as the length of vector C
    sig = sig0 * ones([1,n]);       % sig is defined as a vector of length n with all elements with value sig0
    V = Vega(K, S0, r , q, t, sig); % V is defined as a vector with the vega values for each K and sig value under the Black Scholes model
    f = tol;                        % f is defined as the value of tol
    
    % the Newton method is applied:
    for i = 1:maxIterations
        if(abs(V) < tol)                                        % if all absolute values of V are smaller than tol the loop terminates
            return
        else
            sig = sig - (f ./ (-V));                                        % the vector f ./ (-V) is subtracted to vector sig
            V = Vega(K, S0, r, q, t, sig);                                  % the vector V of Vega values is defined
            [C1, P1, D1, G1, R1, T1] = BlackScholes(S0, K, r, q, t, sig);   % C1 is defined as the Black Scholes estimated Call Option value under the K, S0, r, q, t and sig parameters
            f = C - C1;                                                     % f is defined as the difference between the Call Option price C and C1
            if(abs(f) < tol)
                return
            end
        end
    end
end

