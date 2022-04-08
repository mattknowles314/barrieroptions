%Stock and Options Data
S0 = 100; %Initial stock price
K = 90; %Strike price
B = 130; %Barrier price
r = 0.03; %Risk-free rate
q = 0.05; %Divident yirld
T = 0.5; %Time to maturity (years)

%Simulation Data
Nmc = 100; %Number of simulations

P = 100;
expPrices = zeros(P+1,1);
impPrices = zeros(P+1,1);
crnPrices = zeros(P+1,1);

N_start = 500;
M_start = 1000;

for p=1:P
    expPrices(p) = explicit(S0,K,B,r,q,T,(N+(10*p)),(M+(10*p)));
    impPrices(p) = implicit(S0,K,B,r,q,T,(N+(10*p)),(M+(10*p)));
    crnPrices(p) = crank(S0,K,B,r,q,T,(N+(10*p)),(M+(10*p)));
end

Q = linspace(1,10*P,P-1);

plot(Q,expPrices(1:P-1))
hold on
plot(Q,impPrices(1:P-1))
hold on
plot(Q,crnPrices(1:P-1))
legend('Explicit Scheme','Implicit Scheme','Crank Nicholson')
xlabel("Grid Size increase")
ylabel("Option Price")
hold off