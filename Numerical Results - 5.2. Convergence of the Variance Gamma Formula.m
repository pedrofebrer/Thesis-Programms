%%  NUMERICAL RESULTS - 5.2. CONERGENCE OF THE VARIANCE GAMMA FORMULA  %%


% Variance Gamma Parameters:
C = 1.3574;
G = 5.8704;
M = 14.2699;

% General Parameters 
S0 = 1124.47;
r = 0.019;
q = 0.012;

% Time to Maturity:
tau= [(30 - 18)/365, (7-6)/12+ (30 - 18)/365, (7-3)/12+ (30 - 18)/365, (7)/12+ (30 - 18)/365, (7+3)/12+ (30 - 18)/365, (7+6)/12+ (30 - 18)/365, (7+12)/12+ (30 - 18)/365];

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

% Parameters for the Variance Gamma Convergence:
n=30;

% Double sum series computaion computation:
[CEC, CAC, CC, COPT] = ConvergenceVGSum(C, G, M, S, K, r, q, tau, n);

% Values of the double sum series when n = 10 and 0 <= m , k <= 10:
fprintf('Values of the double sum series when n = 10 and 0 <= m , k <= 10:')
CC(:,:,11,1)

% Conversion of the three Double sum series, for the euclidean norm:
fprintf('Conversion of the three Double sum series, for the euclidean norm:')
CEC