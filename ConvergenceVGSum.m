function [CEC, CAC, CC, COPT] = ConvergenceVGSum(C, G, M, S, K, r ,q, tau,  n1)
    
    % Inputs:
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % S - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % tau - time to maturity vector
    % n1 - number of sums on both variables of the double sum series
    
    % Outputs:
    % COPT - European Call Option formula term values matrix
    % CC -  double sum series value matrix, for each K and tau element, when n is constant for CC(:,:,:,2), when m is constant for CC(:,:,:,1) and when k is constant for CC(:,:,:,3)
    % CEC - the the matrix with the values of the three double sum series after applying the euclidean metric for all 189 cases
    % CAC - the matrix with the values of the three double sum series after summing all the absolute values for all 189 cases
    
    

    mu = log(CharacteristicFunctionVG(-i, C, G, M, 1)); % mu is defined as the log(E[e^(tau X_VG)]),
                                                        % i.e. as the real asset rate of return under 'P' driven by the Variance Gamma process
    
    t1 = length(tau);                       % t1 is defined as the length of vector tau
    r1 = length(K);                         % r1 is defined as the length of vector K
    lo = log(S./K) + (r-q).*tau - mu.*tau;  % lo is defined as the matrix of the values of variable [log] under the given parameters
    Mat = zeros(r1,t1,n1+1,n1+1,n1+1);      % Mat is defined as zero five dimensional matrix
    
    % The values of the terms in C^1_VG for each strike price (in vector K)
    % each time to maturnity (in vector tau) and values between 1 and n1, 1
    % and m,1 and 1 and k1, are computed in the following 5 "for" loops
    
    for t = 1:t1
        for l = 1:r1
            for n = 0:n1+1
                for m = 0:n1+1
                    for k = 0:n1+1
                        Mat(l,t,n+1,m+1,k+1) = ((-1)^(n+m)) * gamma(C*tau(t) + m) * gamma(-1 - 2*C*tau(t) -k-n-m) * (M^n) * (G^m) * (- lo(l,t) )^(1 + 2*C*tau(t) +k+n+m) / (gamma(1 - C*tau(t) -n)*factorial(n)*factorial(m));  % Mat(l,t,n+1,m+1,k+1) is defined as the value of the term in C^1_VG for 
                    end                                                                                                                                                                                                         % the corresponding K(l), tau(t), n, m and k computed according to the formula   
                end
            end
        end
    end

    
    Nat = zeros(r1,t1,n1+1,n1+1,n1+1);  % Nat is defined as zero five dimensional matrix
    
    % The values of the terms in C^2_VG for each strike price (in vector K)
    % each time to maturnity (in vector tau) and values between 1 and n1, 1
    % and m,1 and 1 and k1, are computed in the following 5 "for" loops
    
    for t = 1:t1
        for l = 1:r1
            for n = 0:n1+1
                for m = 0:n1+1
                    for k = 0:n1+1
                        Nat(l,t,n+1,m+1,k+1) = ((-1)^(n+m)) * gamma(C*tau(t) + m) * gamma(1 + 2*C*tau(t) +k-n+m) * (M^(-1-2*C*tau(t)+n-k-m)) * (G^m) * (- lo(l,t))^n / (gamma(2 + C*tau(t) -n+k+m)*factorial(n)*factorial(m));  % Nat(l,t,n+1,m+1,k+1) is defined as the value of the term in C^2_VG for 
                    end                                                                                                                                                                                                         % the corresponding K(l), tau(t), n, m and k computed according to the formula
                end
            end
        end
    end

    Kat = zeros(r1,t1, n1+1, n1+1,n1+1);   % Kat is defined as zero five dimensional matrix
    
    % The values of the terms in C^3_VG for each strike price (in vector K)
    % each time to maturnity (in vector tau) and values between 1 and n1, 1
    % and m,1 and 1 and k1, are computedin the following 5 "for" loops
    
    for t = 1:t1
        for l = 1:r1
            for n = 0:n1+1
                for m = 0:n1+1
                    for k = 0:n1+1
                        if  (lo(l,t) > 0)
                            Kat(l,t,n+1,m+1,k+1) = ((-1)^(C*tau(t) + m)) * gamma(C*tau(t) + n) * gamma(C*tau(t) + m) * (M^n) * (G^m) * ((lo(l,t))^(1 + 2*C*tau(t) +k+m+n)) / (gamma(2 + 2*C*tau(t) +k+n+m)*factorial(n)*factorial(m));  % if [log] > 0 Kat(l,t,n+1,m+1,k+1) is defined as the value of the term in C^3_VG for 
                        else                                                                                                                                                                                                            % the corresponding K(l), tau(t), n, m and k computed according to the formula,
                            Kat(l,t,n+1,m+1,k+1) = 0;                                                                                                                                                                                   % else if [log] < 0 Kat(l,t,n+1,m+1,k+1) is defined as 0
                        end
                    end                                                                                                                                                                                                                  
                end                                                                                                                                                                                                                     
            end
        end
    end 
    
    COPT =    real( K .* ((G*M).^(C*tau)) .* exp( -(r-q)*tau) .* (Mat + Nat + Kat) ./ ((gamma(C*tau)).^2));  % COPT is defined as the five dimensional matrix for the terms of the triple sum
                                                                                                            % series formula for an European Call Option formula for each element of K and tau
    
    CC = zeros(r1,t1,n1+1,3);   % CC is defined as a three dimensional zero matrix
    
    for n = 1:n1+1
        for j = 1:n
            for k = 1:n
                CC(:,:,n,1) = CC(:,:,n,1) + COPT(:,:,n,j,k); % CC(:,:,:,1) is defined as the double sum series value matrix, for each K and tau element, when n is constant
                CC(:,:,n,2) = CC(:,:,n,2) + COPT(:,:,j,n,k); % CC(:,:,:,2) is defined as the double sum series value matrix, for each K and tau element, when m is constant
                CC(:,:,n,3) = CC(:,:,n,3) + COPT(:,:,j,k,n); % CC(:,:,:,3) is defined as the double sum series value matrix, for each K and tau element, when k is constant
            end
        end
    end
    
    CeC =permute(sqrt(sum(sum(CC.^2,2),1)),[3 4 1 2]);
    CEC = vpa([0:n1 ; CeC']',8);                                % CEC is defined as the matrix with the values of the double sum series after applying the euclidean metric for all 189 cases
    CaC = permute(sum(sum(abs(CC),2),1),[3 4 1 2]);             
    CAC = vpa([0:n1 ; CaC']',8);                                % CAC is defined as the matrix with the values of the double sum series after summing all the absolute values for all 189 cases
end