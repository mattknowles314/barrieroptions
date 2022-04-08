%Stock and Options Data
S0 = 100; %Initial stock price
K = 90; %Strike price
B = 130; %Barrier price
r = 0.03; %Risk-free rate
q = 0.05; %Divident yirld
T = 0.5; %Time to maturity (years)

%Simulation Data
Nmc = 100; %Number of simulations
M = 160; %Number of time steps
N = 150; %Number of space steps

[MCPrice, MCError] = monte_carlo(S0,K,B,r,q,T,Nmc,M);
EXPPrice = explicit(S0,K,B,r,q,T,N,M);
IMPPrice = implicit(S0,K,B,r,q,T,N,M);

MCPrice
EXPPrice
IMPPrice