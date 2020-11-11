%%  NUMERICAL RESULTS - 5.4. THE GREEKS FORMULAS BEHAVIOR  %%


%% DELTA VARIANCE GAMMA


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters:
S = 500:5:2000;
K = 1025;
tau=35*7/365;
r = 0.019;
q = 0.012;

% Stochastic Volatility:
sig = 0.1812;

% Parameters for the Variance Gamma Formula:
n=32;
m=19;
k=6;


% Delta measure for the Black-Scholes model computation:
[CallBS, PutBS, DeltaBS, GammaBS, RhoBS, ThetaBS] = BlackScholes(S, K, r, q, tau, sig);

% Delta measure for the Variance Gamma formula computation:
len = length(S);
DD = zeros(len,1);
for j = 1:len
    Delta = DeltaVGSumFormula(C, G, M, S(j), K, r ,q, tau,  n, m, k);
    DD(j) = Delta;
end


% Plot of the Delta measure function
figure;
hold on
grid off
plot(S,DeltaBS,'r',S,DD,'b','LineWidth',0.1);
xlabel 'Stock Price';
ylabel 'Variance gamma values';
title 'Delta';
legend({'Brownian Motion','Variance Gamma Formula'}, 'Location', 'NorthWest');
hold off

%% GAMMA VARIANCE GAMMA


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters:
S = 500:5:2000;
K = 1025;
tau=35*7/365;
r = 0.019;
q = 0.012;

% Volatility:
sig = 0.1812;

% Parameters for the Variance Gamma Formula:
n=38;
m=27;
k=7;


% Gamma measure for the Black-Scholes model computation:
[CallBS, PutBS, DeltaBS, GammaBS, RhoBS, ThetaBS] = BlackScholes(S, K, r, q, tau, sig);

% Gamma measure for the Variance Gamma formula computation:
len = length(S);
GG = zeros(len,1);
for j = 1:len
    Gamma = GammaVGSumFormula(C, G, M, S(j), K, r ,q, tau,  n, m,k);  
    GG(j) = Gamma;
end


% Plot of the Gamma measure function
figure;
hold on
grid off
plot(S,GammaBS,'r',S,GG,'b','LineWidth',0.1);
xlabel 'Stock Price (USD)';
ylabel 'Gamma values';
title 'Gamma';
legend({'Brownian Motion','Variance Gamma Formula'}, 'Location', 'NorthEast');
hold off

%% RHO VARIANCE GAMMA


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters:
S = 500:5:2000;
K = 1025;
tau=35*7/365;
r = 0.019;
q = 0.012;

% Volatility:
sig = 0.1812;

% Parameters for the Variance Gamma Formula:
n=38;
m=27;
k=8;


% Rho measure for the Black-Scholes model computation:
[CallBS, PutBS, DeltaBS, GammaBS, RhoBS, ThetaBS] = BlackScholes(S, K, r, q, tau, sig);

% Rho measure for the Variance Gamma formula computation:
len = length(S);
RR = zeros(len,1);
for j = 1:len
    Rho = RhoVGSumFormula(C, G, M, S(j), K, r ,q, tau,  n, m,k);
    RR(j) = Rho;
end


% Plot of the Rho measure function
figure;
hold on
grid off
plot(S,RhoBS,'r',S,RR,'b','LineWidth',0.1);
xlabel 'Stock Price (USD)';
ylabel 'Rho values';
title 'Rho';
legend({'Brownian Motion','Variance Gamma Formula'},'Location', 'NorthWest');
hold off

%% THETA VARIANCE GAMMA


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters
S = 500:5:1500;
K = 1025;
tau=35*7/365;
r = 0.019;
q = 0.012;

% Volatility:
sig = 0.1812;

% Parameters for the Variance Gamma formula:
n=38;
m=27;
k=7;

% Parameter for the Digamma function
p=20000;


% Theta measure for the Black-Scholes model computation:
[CallBS, PutBS, DeltaBS, GammaBS, RhoBS, ThetaBS] = BlackScholes(S, K, r, q, tau, sig);

% Theta measure for the Variance Gamma formula computation:
len = length(S);
TT = zeros(len,1);
for j = 1:len
    S(j)
    Theta = ThetaVGSumFormula(C, G, M, S(j), K, r ,q, tau,  n, m, k, p);
    TT(j) = Theta;
end


% Plot of the Theta measure function
figure;
hold on
grid off
%plot(S,ThetaBS,'r','LineWidth',0.1);
plot(S,ThetaBS,'r',S,TT,'b','LineWidth',0.1);
xlabel 'Stock Price (USD)';
ylabel 'Theta Values';
title 'Theta';
legend({'Brownian Motion','Variance Gamma Formula'}, 'Location', 'NorthEast');
hold off