Smin = 0;
Smax = 4*K;

N = 800;
S = linspace(Smin,Smax,N+1);
dS = S(3)-S(2);

T = 0.5;
M = 1000;
tau = linspace(0,T,M+1);
dtau = tau(2)-tau(1);

sigs = zeros(N+1,M);
for j = 2:M
    for k = 2:N
        sigs(j,k) = 0.25*exp(-tau(j))*(100/S(k))^0.35;
    end
end
    
%Matrices for storing values of Alpha and Beta at a given point
Alpha = zeros(M,N+1);
Beta = zeros(M,N+1);
    
%Calculate all values for alpha and beta
for j=1:N+1
    for k=1:M
        Alpha(j,k) = 0.5*(sigs(j,k)^2)*(S(j))^2*(dtau/(dS^2));
        Beta(j,k) = (r-q)*S(j)*(dtau/(2*dS));
    end
end

L = zeros(M+1,N+1);
U = zeros(M+1,N+1);
D = zeros(M+1,N+1);

for j=2:M
    for k=1:N
        L(j,k+1) = -Alpha(j,k+1)+Beta(j,k+1);
        U(j,k+1) = -Alpha(j,k+1)-Beta(j,k+1);
        D(j,k+1) = 1+r*dtau + 2*Alpha(j,k+1);
    end
end

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

%Solve Black-Scholes

for k=1:N
        Vnew(k) = L(1)*V(k+1) + D(j,k+1)*V(j+1) + U(j,k+1)*V(j+1);
end

