function [T1,X,s1] = VarianceGamma(C, G, M, T, dt)
    
    % Inputs:
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % T - time vector
    % dt - time step
    
    % Outputs:
    % T1 - time matrix
    % X - Variance Gamma process simulation matrix for each time value in vector T
    % s1 - vector which contains the number of preformed time steps by the simulation for each time value of in vector T

    [T1,G1,s1] = GammaProcess(C, M, T, dt); % G1 is defined as the Gamma Process simulation matrix for each time value in vector T, for the parameters C, M
    [T2,G2,s2] = GammaProcess(C, G, T, dt); % G2 is defined as the Gamma Process simulation matrix for each time value in vector T, for the parameters C, G
    
    X = G1 - G2;    % X is defined as the difference between the two Gamma process simulation matrices G1 and G2
                    % i.e. as th Variance Gamma process simulation matrix for each time value in vector T, for the parameters G, G, M
    
end