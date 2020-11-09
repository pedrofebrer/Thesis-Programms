function phi = CharacteristicFunctionVG(u, C, G, M, T)
    
    % Inputs:
    % u - characteristic function input variable
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % T - time to maturity vector

    
    % Outputs:
    % phi - value of the characteristic function
    
    phi = (G.*M ./ (G.*M + (M - G).*i.*u + u.^2)).^(C.*T);  % phi is defined by the value for the formula for the
                                                            % characteristic function of a Variance Gamma process
    
end