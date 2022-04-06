%Stock and Options Data
S0 = 100; %Initial stock price
K = 90; %Strike price
B = 130; %Barrier price
r = 0.03; %Risk-free rate
q = 0.05; %Divident yirld
T = 0.5; %Time to maturity (years)

%Simulation Data
Nmc = 100; %Number of simulations
M = 500; %Number of time steps
N = 500;

%Test Data
%Ntests = 100;
%priceVals = zeros(1,Ntests);
%errorVals = zeros(1,Ntests);
%for n = 1:Ntests
%    [price,error] = monte_carlo(S0,K,B,r,q,T,n,M);
%    priceVals(n) = price;
%    errorVals(n) = error;
%end

%plot(errorVals)
%title("Standard error, with antithetic sampling")
%xlabel("Number of samples in Monte Carlo simmulation")
%ylabel("Standard Error")


%[price] = explicit(S0,K,B,r,q,T,N,M);