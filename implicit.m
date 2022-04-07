function [price] = implicit(S0,K,B,r,q,T,N,M)
    Smin=0;
    Smax=4*K;

    S = linspace(Smin,Smax,N+1);
    dS = S(2) - S(1);
    
    tau = linspace(0,T,M+1);
    dtau = tau(2) - tau(1);
    
    L = ones(M+1,N+1);
    U = ones(M+1,N+1);
    D = ones(M+1,N+1);
    
    diffTau1 = dtau/(dS)^2;
    diffTau2 = dtau/(2*dS);
    
    for j=2:M
        for k=2:N+1
            %Calculate local volatility
            sig = 0.25*exp(-tau(j))*(100/S(k))^0.35;
    
            %Calculate alpha and beta
            alpha = 0.5*(sig^2)*(S(k)^2)*diffTau1;
            beta = (r-q)*S(k)*diffTau2;
            
            %Calculate diagonals
            L(j,k) = -alpha+beta;
            D(j,k) = 1+r*dtau +2*alpha;
            U(j,k) = -alpha-beta;
        end
    end
    
    %Empty array for options prices
    V = zeros(N+1,1); 
    %Initial condition, V(S,0) = max(S-K,0) if S<B, 
    for k = 1:N-1
        if S(k) < B
            V(k) = max(S(k)-K,0);
        else
            V(k) = 0;
        end
    end
    Vnew = V;
    
    %Solve black scholes
    for j=1:N-1
        %Turn RHS into one matrix
        VRHS = Vnew+[-L(2,j+1)*V(1);zeros(N-1,1);-U(M,j+1)*V(N+1)];
        %Solve tridiagonal step
        Vnew = tridiag(D(j+1,:),U(j+1,:),L(j+1,:),VRHS);
    end
    %Update the option prices
    V = Vnew;
    
    price = interp1(S,V,S0);
end

%Tri-diagonal Matrix Solver - Modified from Dr. Francesco Cosentino's
%original algorithm
function x = tridiag (Dx,Ux,Lx,B)
  n = length(B);
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

