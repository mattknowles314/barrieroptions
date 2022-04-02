S0 = 100; %Initial stock price
K = 90; %Strike price
B = 130; %Barrier price
r = 0.03; %Risk-free rate
q = 0.05; %Divident yirld
T = 0.5; %Time to maturity (years)
M=1000;

tic
    [price,error] = monte_carlo(S0,K,B,r,q,T,100000,M);
toc