S0 = 100;
K = 90;
B = 130;
r = 0.03;
q = 0.05;
T = 0.5;

Smin = 0;
Smax = 4*K;

%Price Grid
    
S = linspace(Smin,Smax,N+1);
dS = S(2)-S(1);

%Time Grid
tau = linspace(0,T,M+1);
dtau = tau(2)-tau(1);
    
%Empty array for options prices
V = zeros(N+1,1); 
%Initial condition, V(S,0) = max(S-K,0) if S<B, 
for p = 1:N+1
    if S(p) < B
        V(p) = max(S(p)-K,0);
    else
        V(p) = 0;
    end
end
Vnew = V;

%Tridiagonals for implicit scheme
LI = ones(M+1,N+1);
UI = ones(M+1,N+1);
DI = ones(M+1,N+1);

LE = ones(M+1,N+1);
UE = ones(M+1,N+1);
DE = ones(M+1,N+1);
    
diffTau1 = dtau/(dS)^2;
diffTau2 = dtau/(2*dS);
    
for j=2:M
   for k=2:N+1
       %Calculate local volatility
       sig = 0.25*exp(-tau(j))*(100/S(k))^0.35;
    
       %Calculate alpha and beta
       alpha = 0.5*(sig^2)*(S(k)^2)*diffTau1;
       beta = (r-q)*S(k)*diffTau2;
            
       %Calculate diagonals for implicit scheme
       LI(j,k) = -alpha+beta;
       DI(j,k) = 1+r*dtau+2*alpha;
       UI(j,k) = -alpha-beta;

       %Calculate diagonals for explicit scheme
       LE(j,k) = alpha-beta;
       DE(j,k) = 1-r*dtau-2*alpha;
       UE(j,k) = alpha+beta;
   end


   
end