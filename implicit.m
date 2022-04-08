function [price] = implicit(S0,K,B,r,q,T,N,M)
%Stock Data
    Smin = 0;
    Smax = 4*K;
    
    %Price Grid
    
    S = linspace(Smin,Smax,N+1);
    dS = S(2)-S(1);
    S1 = S(2:N+1);

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
            L(j,k) = -Alpha(j,k)+Beta(j,k);
            D(j,k) = 1+(r*dtau)+(2*Alpha(j,k));
            U(j,k) = -Alpha(j,k)-Beta(j,k);
        end
    end
    
    %Solve Black-Scholes PDE
    for k=1:M-1
        D1 = D(k+1,1) + 2*L(k+1,1);
        DN = D(k+1,N) + 2*U(k+1,N);
        U1 = U(k+1,1) - L(k+1,1);
        LN = L(k+1,N) - U(k+1,N);
        
        Dx = [D1,D(k,1:N-1),DN];
        Ux = [U1,U(k,1:N)];
        Lx = [L(k,1:N),LN];
        
        Vnew = tridiag(Dx,Ux,Lx,V);
    end
    V = Vnew;
    
    %Use 1D Interpolation to price the option at S0
    price = interp1(S1,V,S0);
end

function x = tridiag (Dx,Ux,Lx,B)
  n = length(B)-1;
  x = zeros(n,1);

  %make A upper diagonal
  for p=2:n
      Dx(p) = Dx(p) - Lx(p)*Ux(p-1) / Dx(p-1);
      B(p) = B(p) - Lx(p)*B(p-1) / Dx(p-1);
  end
  
  %Back substitution
  x(n) = B(n) / Dx(n);

  for p=n-1:-1:1
      x(p) = (B(p) - Ux(p)*x(p+1)) / Dx(p);
  end
end

