%%  NUMERICAL RESULTS - 5.3. OPTION PRICE BEHAVIOR FOR THE VARIANCE GAMMA FORMULA  %%


%% EUROPEAN CALL AND PUT OPTION PRICES VS STOCK VALUE


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters:
S = 500:5:4000;
K = 1110;
r = 0.019;
q = 0.012;
tau = 7/12 + (30 - 18)/365;

% Volatility:
sig = 0.1812;

% Parameters for the Variance Gamma Formula:
n=58;
m=31;
k=11;


% Call and Put Option prices for the Black-Scholes model computation:
[CBS, PBS, Delta, Gamma, Rho, Theta] = BlackScholes(S, K, r, q, tau, sig);


% Call and Put Option prices for the Variace Gamma formula computation:
CVGF = zeros(1,length(S));
PVGF = zeros(1,length(S));
for s =1:length(S)
    [COP,POP,lo,mu] = EuropeanOptionSumFormula(C, G, M, S(s), K, r ,q, tau, n, m,k);
    CVGF(s) = COP;
    PVGF(s) = POP;
end


% Plot of the Call Option Price:
figure;
hold on
grid off
plot(S, CVGF,'b',S, CBS,'r--');
xlabel 'Asset Price (USD)';
ylabel 'Option Price (USD)';
title 'European Call Option Price Estimate for the S&P500';
legend({'Variance Gamma Formula','Brownian Motion'}, 'Location', 'NorthWest');
hold off

% Plot of the Put Option Prices:
figure;
hold on
grid off
plot(S, PVGF,'b',S, PBS,'r--');
xlabel 'Asset Price (USD)';
ylabel 'Option Price (USD)';
title 'European Put Option Price Estimate for the S&P500';
legend({'Variance Gamma Formula','Brownian Motion'}, 'Location', 'NorthEast');
hold off

%% VOLATILITY SMILE


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters:
S0 = 1124.47;
tau = 7/12 + (30 - 18)/365;
r = 0.019;
q = 0.012;

%Strike Prices:
K = 975:5:1275;
KData = [975 995 1025 1050 1075 1100 1125 1140 1150 1175 1200 1225 1250 1275];
K1 = K';

% Volatility:
sig = 0.1812;

% Parameters for the Variance Gamma Formula:
n=21;
m=24;
k=7;

% Parameters for the Monte Carlo Simulation:
dt=0.0001;
ite = 10000;

% Actual Values of the Call Option:
Data = [173.30 157.00 133.10 114.80 97.60 81.20 66.90 58.90 53.90 42.50 33.00 24.90 18.30 13.20];

% Call Option prices for the Variace Gamma formula computation:
[CC1,PP1,lo,mu] = EuropeanOptionSumFormula(C, G, M, S0, K1, r ,q, tau, n, m,k);
CVGF = CC1';

% Call Option prices for the Black-Scholes model computation:
[CBS, PBS, Delta, Gamma, Rho, Theta] = BlackScholes(S0, K, r, q, tau, sig);

% Call Option prices for the Variace Gamma Monte Carlo simulation:
CC2 = CallOptionVGRevenge(C, G, M, S0, K, r ,q, tau, dt, ite);
CVGMC = CC2';


% Plot of the Call Option Values:
figure;
hold on;
grid off;
plot(K,CBS, 'r',K,CVGMC,'g', K, CVGF, 'b',KData,Data,'ko','LineWidth',0.1);
xlabel 'Strike Price (USD)';
ylabel 'Value (USD)';
title 'European Call Option Price for the S&P500';
legend({'Black Scholes', 'Variance Gamma Monte Carlo Simulation', 'Variance Gamma with Stochastic Volatility','Observed Prices'}, 'Location', 'NorthEast');
hold off

% Parameters for the Newton method:
sig0 = 0.5;
tol = 0.000001;
maxIterations = 1000;

% Implied Volatilities:
% Implied Volatility for the Black-Scholes computation:
sigBS = ImpliedVolatility (CBS, K, S0, r , q, tau, sig0, tol, maxIterations);
% Implied Volatility for the Variance Gamma Monte Carlo simulation computation:
sigVGMC = ImpliedVolatility (CVGMC, K, S0, r , q, tau, sig0, tol, maxIterations);
% Implied Volatility for the Variance Gamma formula computation:
sigVGF = ImpliedVolatility (CVGF, K, S0, r , q, tau, sig0, tol, maxIterations);
% Implied Volatility for the empirical results:
sigData = ImpliedVolatility (Data, KData, S0, r , q, tau, sig0, tol, maxIterations);


% Plot of the Implied Volatilities:
figure;
hold on
grid off
plot(K, sigBS, 'r', K, sigVGMC,'g', K, sigVGF, 'b',KData,sigData,'ko','LineWidth',0.1);
xlabel 'Strike Price (USD)';
ylabel 'Volatility (USD)';
title 'Implied Volatility for each of the Processes';
legend({'Black Scholes', 'Variance Gamma Monte Carlo Simulation', 'Variance Gamma Formula','Observed Volatility'}, 'Location', 'NorthEast');
hold off