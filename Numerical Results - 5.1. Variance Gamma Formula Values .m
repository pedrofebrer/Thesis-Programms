%%  NUMERICAL RESULTS - 5.1. VARIANCE GAMMA FORMULA VALUES %%


%%  VARIANCE GAMMA FORMULA  %%


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters 
S0 = 1124.47;
r = 0.019;
q = 0.012;

% Times to Maturity:
tau=[4*7/365  9*7/365  22*7/365  35*7/365  48*7/365  61*7/365  87*7/365];

% Strike Prices:
K =[975;
    995;
    1025;
    1050;
    1075;
    1090;
    1100;
    1110;
    1120;
    1125;
    1130;
    1135;
    1140;
    1150;
    1160;
    1170;
    1175;
    1200;
    1225;
    1250;
    1275;
    1300;
    1325;
    1350;
    1400;
    1450;
    1500];

% Parameters for the Variance Gamma formula:
n=22;
m=27;
k=7;

% Call Option prices for the Variace Gamma formula computation:
[COP,POP,lo,mu]  = VarianceGammaSumFormula(C, G, M, S0, K, r ,q, tau, n, m, k);
fprintf('Call Option prices for the Variace Gamma formula:')
COP

% Actual Values of the Call Option:
Data = [     0      0 161.60 173.30      0      0      0;
             0      0 144.80 157.00      0 182.10      0;
             0      0 120.10 133.10 146.50      0      0;
             0  84.50 100.70 114.80      0 143.00 171.40;
             0  64.30  82.50  97.60      0      0      0;
         43.10      0      0      0      0      0      0;
         35.60      0  65.50  81.20  96.20 111.30 140.40;
             0  39.50      0      0      0      0      0;
         22.90  33.50      0      0      0      0      0;
         20.20  30.70  51.00  66.90  81.70  97.00      0;
             0  28.00      0      0      0      0      0;
             0  25.60  45.50      0      0      0      0;
         13.30  23.20      0  58.90      0      0      0;
             0  19.10  38.10  53.90  68.30  83.30 112.80;
             0  15.30      0      0      0      0      0;
             0  12.10      0      0      0      0      0;
             0  10.90  27.70  42.50  56.60      0  99.80;
             0      0  19.60  33.00  46.10  60.90      0;
             0      0  13.20  24.90  36.90  49.80      0;
             0      0      0  18.30  29.30  41.20  66.90;
             0      0      0  13.20  22.50      0      0;
             0      0      0      0  17.20  27.10  49.50;
             0      0      0      0  12.80      0      0;
             0      0      0      0      0  17.10  35.70;
             0      0      0      0      0  10.10  25.20;
             0      0      0      0      0      0  17.00;
             0      0      0      0      0      0  12.20];
             
            
% Root Mean Square Error for the Variace Gamma formula:
err = RMSE(Data,COP);
fprintf('Root Mean Square Error for the Variace Gamma formula:')
err


%%  VARIANCE GAMMA MONTE CARLO SIMULATION  %%


% General Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters 
S0 = 1124.47;
r = 0.019;
q = 0.012;

% Time to Maturity:
tau=[4*7/365  9*7/365  22*7/365  35*7/365  48*7/365  61*7/365  87*7/365];

% Strike Prices:
K =[975;
    995;
    1025;
    1050;
    1075;
    1090;
    1100;
    1110;
    1120;
    1125;
    1130;
    1135;
    1140;
    1150;
    1160;
    1170;
    1175;
    1200;
    1225;
    1250;
    1275;
    1300;
    1325;
    1350;
    1400;
    1450;
    1500];

% Parameters for the Monte Carlo Simulation:
dt=0.0001;
ite = 10000;

% Call Option prices for the Variace Gamma Monte Carlo simulation:
COP = CallOptionMonteCarlo(C, G, M, S0, K, r ,q, tau, dt, ite);
fprintf('Call Option prices for the Variace Gamma Monte Carlo simulation:')
COP

% Actual Values of the Call Option:
Data = [     0      0 161.60 173.30      0      0      0;
             0      0 144.80 157.00      0 182.10      0;
             0      0 120.10 133.10 146.50      0      0;
             0  84.50 100.70 114.80      0 143.00 171.40;
             0  64.30  82.50  97.60      0      0      0;
         43.10      0      0      0      0      0      0;
         35.60      0  65.50  81.20  96.20 111.30 140.40;
             0  39.50      0      0      0      0      0;
         22.90  33.50      0      0      0      0      0;
         20.20  30.70  51.00  66.90  81.70  97.00      0;
             0  28.00      0      0      0      0      0;
             0  25.60  45.50      0      0      0      0;
         13.30  23.20      0  58.90      0      0      0;
             0  19.10  38.10  53.90  68.30  83.30 112.80;
             0  15.30      0      0      0      0      0;
             0  12.10      0      0      0      0      0;
             0  10.90  27.70  42.50  56.60      0  99.80;
             0      0  19.60  33.00  46.10  60.90      0;
             0      0  13.20  24.90  36.90  49.80      0;
             0      0      0  18.30  29.30  41.20  66.90;
             0      0      0  13.20  22.50      0      0;
             0      0      0      0  17.20  27.10  49.50;
             0      0      0      0  12.80      0      0;
             0      0      0      0      0  17.10  35.70;
             0      0      0      0      0  10.10  25.20;
             0      0      0      0      0      0  17.00;
             0      0      0      0      0      0  12.20];
             
            
% Root Mean Square Error for the Variace Gamma Monte Carlo simulation:
err = RMSE(Data,COP);
fprintf('Root Mean Square Error for the Variace Gamma Monte Carlo simulation:')
err

%%  EXAMPLE OF AN OPTION PRICE SIMULATION VALUES  %%


% One of the results for the Call option price under 
% the Variancr Gamma Monte Carlo simulation:
fprintf('Call Option prices for the Monte Carlo simulation example:')
COP = [151.9653  157.3531  167.5363  176.0217  185.6831  191.3072  209.9063;
       132.6797  138.6938  150.3016  160.0207  170.3483  176.5546  196.0470;
       104.0542  111.2052  125.2270  136.9252  148.2588  155.3427  176.1318;
        80.5440   88.9269  105.2361  118.6385  130.7480  138.5410  160.3732;
        57.4830   67.3345   86.2784  101.3093  114.1550  122.6130  145.3532;
        43.9917   54.8873   75.4615   91.4389  104.6675  113.4960  136.7311;
        35.1832   46.8599   68.5158   85.0639   98.5567  107.6130  131.1635;
        26.5715   39.0817   61.7813   78.8586   92.6186  101.8944  125.7277;
        18.2403   31.6025   55.2949   72.8620   86.8539   96.3553  120.4304;
        14.2489   27.9979   52.1581   69.9435   84.0404   93.6445  117.8313;
        10.4082   24.5153   49.0974   67.0804   81.2735   90.9728  115.2676;
         6.9185   21.1685   46.1129   64.2720   78.5591   88.3516  112.7415;
         5.8351   18.0062   43.2164   61.5221   75.8955   85.7797  110.2471;
         4.5594   12.7935   37.7029   56.1971   70.7278   80.7531  105.3551;
         3.6943   10.2125   32.5900   51.1332   65.7604   75.8923  100.6078;
         3.0543    8.4157   27.9128   46.3362   61.0036   71.2050   95.9982;
         2.7922    7.6875   25.7760   44.0523   58.7011   68.9336   93.7474;
         1.8397    5.0932   17.7131   33.7982   48.0240   58.3183   83.0197;
         1.2832    3.5115   12.7309   25.5572   38.7722   48.8643   73.1726;
         0.9398    2.4881    9.4122   19.4176   31.0259   40.5516   64.1935;
         0.7105    1.8385    7.1123   14.8373   24.7489   33.4362   56.0106;
         0.5501    1.3768    5.4285   11.4368   19.7126   27.4185   48.6282;
         0.4358    1.0514    4.1747    8.9280   15.7748   22.4312   42.0480;
         0.3461    0.8122    3.2323    7.0336   12.6556   18.3781   36.2791;
         0.2303    0.4904    2.0408    4.4070    8.2540   12.2687   26.8192;
         0.1659    0.3042    1.3212    2.8060    5.4967    8.2208   19.7298;
         0.1224    0.1919    0.9073    1.8340    3.7365    5.5250   14.4998];
         
% Actual Values of the Call Option:
Data = [     0      0 161.60 173.30      0      0      0;
             0      0 144.80 157.00      0 182.10      0;
             0      0 120.10 133.10 146.50      0      0;
             0  84.50 100.70 114.80      0 143.00 171.40;
             0  64.30  82.50  97.60      0      0      0;
         43.10      0      0      0      0      0      0;
         35.60      0  65.50  81.20  96.20 111.30 140.40;
             0  39.50      0      0      0      0      0;
         22.90  33.50      0      0      0      0      0;
         20.20  30.70  51.00  66.90  81.70  97.00      0;
             0  28.00      0      0      0      0      0;
             0  25.60  45.50      0      0      0      0;
         13.30  23.20      0  58.90      0      0      0;
             0  19.10  38.10  53.90  68.30  83.30 112.80;
             0  15.30      0      0      0      0      0;
             0  12.10      0      0      0      0      0;
             0  10.90  27.70  42.50  56.60      0  99.80;
             0      0  19.60  33.00  46.10  60.90      0;
             0      0  13.20  24.90  36.90  49.80      0;
             0      0      0  18.30  29.30  41.20  66.90;
             0      0      0  13.20  22.50      0      0;
             0      0      0      0  17.20  27.10  49.50;
             0      0      0      0  12.80      0      0;
             0      0      0      0      0  17.10  35.70;
             0      0      0      0      0  10.10  25.20;
             0      0      0      0      0      0  17.00;
             0      0      0      0      0      0  12.20];
             
% Root Mean Square Error for the Monte Carlo simulation example:
err = RMSE(Data,COP);
fprintf('Root Mean Square Error for the Monte Carlo simulation example:')
err

%%  BLACK-SCHOLES FORMULA  %%


% General Parameters:
S0 = 1124.47;
r = 0.019;
q = 0.012;

% Volatility:
sig = 0.1812;

% Times to Maturity:
tau=[4*7/365  9*7/365  22*7/365  35*7/365  48*7/365  61*7/365  87*7/365];

% Strike Prices:
K =[975;
    995;
    1025;
    1050;
    1075;
    1090;
    1100;
    1110;
    1120;
    1125;
    1130;
    1135;
    1140;
    1150;
    1160;
    1170;
    1175;
    1200;
    1225;
    1250;
    1275;
    1300;
    1325;
    1350;
    1400;
    1450;
    1500];

% Call Option prices for the Black-Scholes model computation:
[Call, Put, Delta, Gamma, Rho, Theta] = BlackScholes(S0, K, r, q, tau, sig);
fprintf('Call Option prices for the Black-Scholes model:')
Call


% Actual Values of the Call Option:
Data = [     0      0 161.60 173.30      0      0      0;
             0      0 144.80 157.00      0 182.10      0;
             0      0 120.10 133.10 146.50      0      0;
             0  84.50 100.70 114.80      0 143.00 171.40;
             0  64.30  82.50  97.60      0      0      0;
         43.10      0      0      0      0      0      0;
         35.60      0  65.50  81.20  96.20 111.30 140.40;
             0  39.50      0      0      0      0      0;
         22.90  33.50      0      0      0      0      0;
         20.20  30.70  51.00  66.90  81.70  97.00      0;
             0  28.00      0      0      0      0      0;
             0  25.60  45.50      0      0      0      0;
         13.30  23.20      0  58.90      0      0      0;
             0  19.10  38.10  53.90  68.30  83.30 112.80;
             0  15.30      0      0      0      0      0;
             0  12.10      0      0      0      0      0;
             0  10.90  27.70  42.50  56.60      0  99.80;
             0      0  19.60  33.00  46.10  60.90      0;
             0      0  13.20  24.90  36.90  49.80      0;
             0      0      0  18.30  29.30  41.20  66.90;
             0      0      0  13.20  22.50      0      0;
             0      0      0      0  17.20  27.10  49.50;
             0      0      0      0  12.80      0      0;
             0      0      0      0      0  17.10  35.70;
             0      0      0      0      0  10.10  25.20;
             0      0      0      0      0      0  17.00;
             0      0      0      0      0      0  12.20];
             
% Root Mean Square Error for the Black-Scholes model:
err = RMSE(Data,Call);
fprintf('Root Mean Square Error for the Black-Scholes model:')
err
