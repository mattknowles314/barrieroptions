%Explicit FDM Scheme for Black-Scholes PDE

function [price] = explicit(S0,K,B,r,q,T,N,M)
%Stock Data
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
    
    %Matrix for storing volatility at different points on the grid
    sigs = zeros(N+1,M);
    for j = 2:N
        for k = 2:M+1
            sigs(j,k) = 0.25*exp(-tau(k))*(100/S(j))^0.35;
    
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
    
    %Matrices for storing upper and lower diagonals
    L = zeros(M+1,N+1);
    D = zeros(M+1,N+1);
    U = zeros(M+1,N+1);
    
    %Calculate lower, main and upper diagonal values
    for j=1:M
        for k=1:N
            L(j,k) = Alpha(j,k)-Beta(j,k);
            D(j,k) = 1-(r*dtau)-(2*Alpha(j,k));
            U(j,k) = Alpha(j,k)+Beta(j,k);
        end
    end
    
    %Solve Black-Scholes PDE
    for k=1:M
        AE_1 = diag(D(k,:));
        AE_2 = diag(U(k,1:N),1);
        AE_3 = diag(L(k,1:N),-1);
        AE = AE_1 + AE_2 + AE_3;
        boundary = [L(k,1)*V(1);zeros(N-1,1);U(k,N)*V(N+1)];
        
        Vnew = AE*V+boundary;
        
    end
    V = Vnew;
    %Use 1D Interpolation to price the option at S0
    price = interp1(S,V,S0); 