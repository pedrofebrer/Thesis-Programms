function [TM,G,s] = GammaProcess(a, b, T, dt)
        
    % Inputs:
    % a - parameter a of the Gamma distribution
    % b - parameter b of the Gamma distribution
    % T - time vector
    % dt - time step
    
    % Outputs:
    % TM - time matrix
    % G - Gamma process simulation matrix for each time value in vector T
    % s - vector which contains the number of preformed time steps by the simulation for each time value of in vector T
    
    l2 = length(T);         % l2 is defined as the length of vector T  
    l1 = ceil(T(end)/dt);   % l1 is defined as number of of time steps the simulation will preform
    TM = zeros(l1,l2);      % TM is defined as a two dimensional zero matrix
    
    s = zeros(1,l2);            % s is defined as a zero vector of length l2
    for k=1:l2
        s(k) = ceil(T(k)/dt);   % s is defined as a vector which each cell  contains
    end                         % the number of steps that will be preformed by the
                                % simulation for each time value of in vector T
    
    for j= 1:l2
        v0 = 0:dt:T(j);
        v0(end+1) = T(j);       % v0 is defined as vector containing the times after each step of the simulation until time T(j)
        v = v0';                % v is defined as transpose of v0
        TM(1:length(v),j) = v;  % TM is defined as the matrix of containing all the vectors v0 for each element in T
    end                        
    
    DT = zeros(l1,l2);                  % DT is defined as a two dimensional matrix
    DT(1,:) = TM(1,:);                  % with the time steps for each element of TM
    for i = 2:l1
        DT(i,:) = TM(i,:) - TM(i-1,:);
    end
    
    V = a*1 + zeros(l1,l2);     % V is defined as a two dimensional matrix where all elements are a
    
    g = JohnksGammaGenerator(V.*DT,TM) / b;     % g is defined as the matrix with element The Johnks Gamma Generator
                                                       % generated for each element of the matrices V.*DT and TM divided by b
                                                       
    G=[zeros(1,l2) ; cumsum(g)];    % G is defined as the Gamma process simulation for each time value in vector T

end