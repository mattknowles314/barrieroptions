%Explicit FDM Scheme for Black-Scholes PDE

%Stock Data
S0 = 100; %Initial stock price
K = 90; %Strike price
B = 130; %Barrier price
r = 0.03; %Risk-free rate
q = 0.05; %Divident yirld
T = 0.5; %Time to maturity (years)
alpha = 0.35; %Used for calculating sigma
Smin = 0;
Smax = 4*K;

%Price Grid
N = 1000;
S1 = linspace(Smin,Smax,N+1);
dS = S1(2)-S1(1);
S = S1(2:N); %Leave first and last out for boundary conditions

%Time Grid
M = 1000;
tau = linspace(0,T,M+1);
dtau = tau(2) -tau(1);


%Solve Black-Scholes


