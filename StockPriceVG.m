function [X,TM,S,s,m] = StockPriceVG(C, G, M, S0, r ,q, T, dt)
    
        
    % Inputs:
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % S0 - initial stock price value
    % r - risk free rate of return
    % q - dividend rate on r
    % T - time to maturity vector
    % dt - time step
     
    % Outputs:
    % TM - time matrix
    % X - Variance Gamma process simulation matrix for each time value in vector T
    % s - vector which contains the number of preformed time steps by the simulation for each time value of in vector T
    % m - matrix of characteristic function values under a Variance Gamma process
    % S - Stock price simulation matrix for each time value of in vector T

    [TM,X,s] = VarianceGamma(C, G, M, T, dt);       % X is defined as the matrix of Variance Gamma process simulations matrix for each time value in vector T
    m = CharacteristicFunctionVG(-i, C, G, M, TM);  % m is defined as the matrix of characteristic function values under a Variance Gamma process

    S = S0.*exp((r-q)*TM+ X) ./ m;  % S is defined as the stock prince value matrix under measure 'Q'
end