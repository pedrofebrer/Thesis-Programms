function COP = CallOptionMonteCarlo(C, G, M, S0, K, r ,q, T, dt, ite)
    
    % Inputs:
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % S0 - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % T - time to maturity vector
    % dt - time step
    % ite - number of iterations
    
    % Outputs:
    % COP - price of the European Call Option matrix
    % s - vector which contains the number of preformed time steps by the simulation for each time value of in vector T
    % POP - prices of the European Put Option matrix
    % lo - values of the parameter [log] matrix
    % mu - real asset rate of return under 'P' driven by the Variance Gamma process
    
    V = zeros(ite,length(T));               % V is defined as two dimensional zero matrix 
    MM = zeros(length(K),length(T),ite);    % MM is defined as a three dimensional zero matrix
    
    for i = 1:ite
        [X,TM,S,s,m] = StockPriceVG(C, G, M, S0, r ,q, T, dt);   % S is defined as the stock price simulation matrix for each time value of in vector T
                                                                        % s is defined as the vector which contains the number of preformed time steps by the simulation for each time value of in vector T

        for j =1:length(T)
            V(i,j) = S(s(j),j); % V is defined as a two dimension matrix with the last stock value of S
        end                     % for each time value of vector T and each iteration ite
        
        MM(:,:,i) = (V(i,:) - K);   % MM is defined as a two dimension matrix where
    end                             % to every element in S is subtracted the value K


    MM(MM<0) = 0;                   % All negative elements of MM are defined as 0
    
    COP = exp(-r*T).*mean(MM,3);    % COP is defined as mean mean matrix of each iteration ite
                                    %  in MM with the interest rate r deduction
                                    % i.e.: COP = exp(-r*T).*E[(MM-K)^+]
end