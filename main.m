%Stock and Options Data
S0 = 100; %Initial stock price
K = 90; %Strike price
B = 130; %Barrier price
r = 0.03; %Risk-free rate
q = 0.05; %Dividend yirld

%Simulation Data
Nmc = 100; %Number of simulations
T = 0.5; %Time to maturity (years)
M = 140; %Number of time steps
N = 150; %Number of space steps

%Calculate Prices
[MCPrice, MCError] = monte_carlo(S0,K,B,r,q,T,Nmc,M); %Calculates price and standard error by monte carlo method
EXPPrice = explicit(S0,K,B,r,q,T,N,M); %Calculate price by explicit method
IMPPrice = implicit(S0,K,B,r,q,T,N,M); %Calculate price by implicit method
CRNPrice = crank(S0,K,B,r,q,T,N,M); %Calculate price by crank-nicholson method

%Ouptut Prices 
MCPrice 
EXPPrice
IMPPrice
CRNPrice