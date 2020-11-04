function Theta = ThetaVGSumFormula(C, G, M, S, K, r ,q, tau,  n1, m1, k1, p)


    % Inputs:
    % C - parameter C of the Variance Gamma process
    % G - parameter G of the Variance Gamma process
    % M - parameter M of the Variance Gamma process
    % S - initial stock price value
    % K - strike price vector
    % r - risk free rate of return
    % q - dividend rate on r
    % tau - time to maturity vector
    % n1 - number of sums on n variable of the triple sum series for the European Call Option Theta measure driven by the Variance Gamma
    % m1 - number of sums on m variable of the triple sum series for the European Call Option Theta measure driven by the Variance Gamma
    % k1 - number of sums on k variable of the triple sum series for the European Call Option Theta measure driven by the Variance Gamma
    % p - number of sums in the diagamma representation series
    
    % Outputs:
    % Theta - Theta measure matrix of the European Call Option driven by the Variance Gamma process

    mu = log(CharacteristicFunctionVG(-i, C, G, M, 1)); % mu is defined as the log(E[e^(tau X_VG)]),
                                                        % i.e. as the real asset rate of return under 'P' driven by the Variance Gamma process
    
    t1 = length(tau);                           % t1 is defined as the length of vector tau
    r1 = length(K);                             % r1 is defined as the length of vector K
    lo = log(S./K) + (r-q).*tau - mu.*tau;      % lo is defined as the matrix of the values of variable [log] under the given parameters
    m0 = C*log(G*M)-(r-q)-2*C*Digamma(C*tau,p); % m0 is defined as the vector of the first three terms of theta_1, theta_2 and theta_3 of the theta measure of the Variance Gamma formula
    MatTheta = zeros(r1,t1,n1+1,m1+1,k1+1);     % MatTheta is defined as zero five dimensional matrix
    
    % The values of the terms in the first sum of the Theta measure for each strike
    % price (in vector K) each time to maturnity (in vector tau) and values between 1
    % and n1, 1 and m and, 1 and k1, are computed in the following 5 "for" loops
    
    for t = 1:t1
        for l = 1:r1
            for n = 0:n1
                for m = 0:m1
                    for k = 0:k1
                        Mat = ((-1)^(n+m)) * gamma(C*tau(t) + m) * gamma(-1 - 2*C*tau(t) -k-n-m) * (M^n) * (G^m) * (- lo(l,t) )^(1 + 2*C*tau(t) +k+n+m) / (gamma(1 - C*tau(t) -n)*factorial(n)*factorial(m));   % Mat is defined as the value of the term in C^1_VG for the corresponding K(l), tau(t), n, m and k computed according to the Variance Gamma formula                                         
                        m11 = C*Digamma(C*tau+m,p)-2*C*Digamma(-1-2*C*tau-k-n-m,p)+C*Digamma(1-C*tau-n,p);                                                                                                      % m11 is defined as the vector of the second three terms of theta_1 of the theta measure of the Variance Gamma formula
                        m12 = 2*C*log(-lo(l,t))+(1+2*C*tau+k+n+m)*(r-q-mu)/lo(l,t);                                                                                                                             % m12 is defined as the vector of the last term of theta_1 plus a term that accounts for the extra term (r-q-mu)*Delta_1 outside of the expression theta_1 * C^1_VG
                        MatTheta(l,t,n+1,m+1,k+1) = (m0 + m11 + m12)*Mat;                                                                                                                                       % MatTheta(l,t,n+1,m+1,k+1) is defined as the value of the term in the first sum of the Theta measure for the corresponding K(l), tau(t), n, m and k computed according to the Variance Gamma formula                               
                    end
                end
            end
        end
    end


    NatTheta = zeros(r1,t1,n1+1,m1+1,k1+1); % NatTheta is defined as zero five dimensional matrix
    
    % The values of the terms in the second sum of the Theta measure for each strike
    % price (in vector K) each time to maturnity (in vector tau) and values between 1
    % and n1, 1 and m and, 1 and k1, are computed in the following 5 "for" loops    
    
    for n = 0:n1
        for m = 0:m1
            for  k = 0:k1
                for t = 1:t1
                    for l = 1:r1
                        Nat = ((-1)^(n+m)) * gamma(C*tau(t) + m) * gamma(1 + 2*C*tau(t) +k-n+m) * (M^(-1-2*C*tau(t)+n-k-m)) * (G^m) * (- lo(l,t))^n / (gamma(2 + C*tau(t) -n+k+m)*factorial(n)*factorial(m));   % Nat is defined as the value of the term in C^2_VG for the corresponding K(l), tau(t), n, m and k computed according to the Variance Gamma formula 
                        m21 = C*Digamma(C*tau+m,p)+2*C*Digamma(1+2*C*tau+k-n+m,p)-C*Digamma(2+C*tau+k-n+m,p);                                                                                                   % m21 is defined as the vector of the second three terms of theta_2 of the theta measure of the Variance Gamma formula
                        m22 = -2*C*log(M) + n*(r-q-mu)/lo(l,t);                                                                                                                                                 % m22 is defined as the vector of the last term of theta_2 plus a term that accounts for the extra term (r-q-mu)*Delta_2 outside of the expression theta_2 * C^2_VG
                        NatTheta(l,t,n+1,m+1,k+1)= (m0 + m21 + m22)*Nat;                                                                                                                                        % NatTheta(l,t,n+1,m+1,k+1) is defined as the value of the term in the second sum of the Theta measure for the corresponding K(l), tau(t), n, m and k computed according to the Variance Gamma formula 
                    end
                end
            end
        end
    end

    KatTheta = zeros(r1,t1,n1+1,m1+1,k1+1); % KatTheta is defined as zero five dimensional matrix
    
    % The values of the terms in the third sum of the Theta measure for each strike
    % price (in vector K) each time to maturnity (in vector tau) and values between 1
    % and n1, 1 and m and, 1 and k1, are computed in the following 5 "for" loops 

    for n = 0:n1
        for m = 0:m1
            for  k = 0:k1
                for t = 1:t1
                    for l = 1:r1
                        if (lo(l,t) > 0)                                                                                                                                                                                % if [log] > 0 
                            Kat = ((-1)^(C*tau(t) + m)) * gamma(C*tau(t) + n) * gamma(C*tau(t) + m) * (M^n) * (G^m) * ((lo(l,t))^(1 + 2*C*tau(t) +k+m+n)) / (gamma(2 + 2*C*tau(t) +k+n+m)*factorial(n)*factorial(m));   % Kat is defined as the value of the term in C^3_VG for the corresponding K(l), tau(t), n, m and k computed according to the Variance Gamma formula
                            m31 = C*Digamma(C*tau+n,p)+C*Digamma(C*tau+m,p)-2*C*Digamma(2+2*C*tau+k+n+m,p);                                                                                                             % m31 is defined as the vector of the second three terms of theta_3 of the theta measure of the Variance Gamma formula
                            m32 = C*pi*i+2*C*log(lo(l,t))+(1+2*C*tau+k+n+m)*(r-q-mu)/lo(l,t);                                                                                                                           % m32 is defined as the vector of the last two terms of theta_3 plus a term that accounts for the extra term (r-q-mu)*Delta_3 outside of the expression theta_3 * C^3_VG
                            KatTheta(l,t,n+1,m+1,k+1) = (m0 + m31 + m32)*Kat;                                                                                                                                           % KatTheta(l,t,n+1,m+1,k+1) is defined as the value of the term in the third sum of the Theta measure for the corresponding K(l), tau(t), n, m and k computed according to the Variance Gamma formula 
                        else                                                                                                                                                                                            % else if [log] < 0 KatTheta(l,t,n+1,m+1,k+1) is defined as 0              
                            KatTheta(l,t,n+1,m+1,k+1) = 0;
                            %Kat2(l,t,n+1,m+1,k+1) = 0;
                        end
                    end
                end
            end
        end
    end
    
    % Tha matrix for the value of the Theta measure is computed:
    Theta =    real( - K .* ((G*M).^(C*tau)) .* exp( -(r-q)*tau) .* (sum(sum(sum(MatTheta,5),4),3) + sum(sum(sum(NatTheta,5),4),3) + sum(sum(sum(KatTheta,5),4),3) ) ./ ((gamma(C*tau)).^2));   % Theta is defined as a two dimensional matrix, that varies as K and tau vary, by summing
                                                                                                                                                                                                % the terms in MatTheta, NatTheta and KatTheta from 1 to n1, 1 to m1 and 1 to k1, respectively
    
end