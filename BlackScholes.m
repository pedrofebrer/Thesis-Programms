function [Call, Put, Delta, Gamma, Rho, Theta] = BlackScholes(S0, K, r, q, t, sig)
    
    % Inputs:
    % S0 - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % t - time to maturity vector
    % sig - volatility value
    
    % Outputs:
    % Call - price of the European Call Option matrix
    % Put - prices of the European Put Option matrix
    % Delta - values of the Delta measure matrix
    % Gamma - values of the gamma measure matrix
    % Rho - values of the Rho measure matrix
    % Theta - values of the Theta measure matrix
    
    d1 = (log(S0./K) + (r-q+(sig.^2)/2)*t) ./(sig.*sqrt(t));    % d1 is defined as a matrix for the values of d_+ of the
                                                                % Black-Scholes formula for the respective K and t values
                  
    d2 = (log(S0./K) + (r-q-(sig.^2)/2)*t) ./(sig.*sqrt(t));    % d2 is defined as a matrix for the values of d_- of the
                                                                % Black-Scholes formula for the respective K and t values
    
    Call = S0*exp(-q*t).*normcdf(d1) - K.*exp(-r*t).*normcdf(d2);                   % Call is defined as a two dimensional matrix with the values for the Black-Scholes
                                                                                    % formula for an European Call Option for the respective K and t values
    
    Put = K.*exp(-r*t).*normcdf(-d2) - S0.*exp(-q*t).*normcdf(-d1);                 % Put is defined as a two dimensional matrix with the values for the Black-Scholes
                                                                                    % formula for an European Put Option for the respective K and t values
    
    Delta = normcdf(d1);                                                            % Delta is defined as a two dimensional matrix with the values for the Black-Scholes formula
                                                                                    % for the Delta measurement for an European Csll Option for the respective K and t values
    
    Gamma = normpdf(d1)./(S0.*sig.*sqrt(t));                                        % Gamma is defined as a two dimensional matrix with the values for the Black-Scholes formula
                                                                                    % for the Gamma measurement for an European Csll Option for the respective K and t values
    
    Rho = K.*t.*exp(-r.*t).*normcdf(d2);                                            % Rho is defined as a two dimensional matrix with the values for the Black-Scholes formula
                                                                                    % for the Rho measurement for an European Csll Option for the respective K and t values
    
    Theta = - (S0.*normpdf(d1).*sig)./(2*sqrt(t)) - r.*K.*exp(-r.*t).*normcdf(d2);  % Theta is defined as a two dimensional matrix with the values for the Black-Scholes formula
                                                                                    % for the Theta measurement for an European Csll Option for the respective K and t values
    
end