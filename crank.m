function [price] = crank(S0,K,B,r,q,T,N,M)
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
    
    %Solve Black-Scholes PDE
    for k=1:M-1
        %Calculate AE matrix by calculating all diagonals
        AE_1 = diag(DE(k,:),0);
        AE_2 = diag(UE(k,1:N),1);
        AE_3 = diag(LE(k,1:N),-1);    
        AE = AE_1+AE_2+AE_3; %Sum diagonals into one matrix
        
        R = (AE+eye(N+1))*V; %Calculate (AE+I)V term
        boundary = [LE(k,1)*V(1)-LI(k+1,2)*Vnew(N+1);zeros(N-1,1);UE(k,N)*V(N+1)-UI(k+1,N)*Vnew(N+1)];
        VRHS = R + boundary; %This gives the (N+1)x1 vector on the RHS of the CN scheme equation 
    
        AI_D = DI(k,:)+1; %Add one to this term as we have (AI+I)
        AI_U = UI(k,:);
        AI_L = LI(k,:);
    
        Vnew = tridiag(AI_D,AI_U,AI_L,VRHS); %Why is this shrinking?
    end
    V = Vnew;
    
    price = interp1(S,V,S0)
end

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