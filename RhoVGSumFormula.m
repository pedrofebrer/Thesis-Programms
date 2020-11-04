function Rho = RhoVGSumFormula(C, G, M, S, K, r ,q, tau, n1, m1, k1)
    % Inputs:
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % S - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % tau - time to maturity vector
    % n1 - number of sums on n variable of the triple sum series for the European Call Option Rho measure driven by the Variance Gamma
    % m1 - number of sums on m variable of the triple sum series for the European Call Option Rho measure driven by the Variance Gamma
    % k1 - number of sums on k variable of the triple sum series for the European Call Option Rho measure driven by the Variance Gamma
    
    % Outputs:
    % Rho - Rho measure matrix of the European Call Option driven by the Variance Gamma process
    
    D = DeltaVGSumFormula(C, G, M, S, K, r ,q, tau,  n1, m1,k1);                                % D is defined as the two dimensional matrix of the delta measure values for an European Call option driven by the Variance Gamma
    [COP,POP,lo,mu] = EuropeanOptionSumFormula(C, G, M, S, K, r ,q, tau,  n1, m1,k1);           % COP is defined as the two dimensional matrix of the European Call option price driven by the Variance Gamma
    Rho = -tau*COP + tau*S*D;                                                                   % Rho is defined as the two dimensional matrix of the rho measure values for an European Call option driven by the Variance Gamma
end