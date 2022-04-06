%Stock Data
Smin = 0;
Smax = 4*K;
    
%Price Grid
    
S = linspace(Smin,Smax,N+1);
dS = S(2)-S(1);

%Time Grid
tau = linspace(0,T,M+1);
dtau = tau(2)-tau(1);
    
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
    for j = 2:M
        for k = 2:N+1
            sigs(j,k) = 0.25*exp(-tau(j))*(100/S(k))^0.35;
        end
    end
    
    %Matrices for storing values of Alpha and Beta at a given point
    Alpha = zeros(M,N+1);
    Beta = zeros(M,N+1);
    
    %Calculate all values for alpha and beta
    for j=1:N+1
        for k=1:M
            Alpha(j,k+1) = 0.5*(sigs(j,k+1)^2)*(S(j))^2*(dtau/(dS^2));
            Beta(j,k+1) = (r-q)*S(j)*(dtau/(2*dS));
        end
    end
    
    %Matrices for storing upper and lower diagonals
    L = zeros(M,N+1);
    D = zeros(M,N+1);
    U = zeros(M,N+1);

    %Calculate lower, main and upper diagonal values
    for j=1:M
        for k=1:N
            L(j,k) = -Alpha(j,k+1)+Beta(j,k+1);
            D(j,k) = 1+(r*dtau)+(2*Alpha(j,k+1));
            U(j,k) = -Alpha(j,k+1)-Beta(j,k+1);
        end
    end 
    
    %These values are used in the diagomal matrix later on
    Dprime1 = zeros(N+1,1);
    Uprime1 = zeros(N+1,1);
    DprimeN = zeros(N+1,1);
    LprimeN = zeros(N+1,1);

    for k=1:N
        Dprime1(k+1) = D(1,k+1)+2*L(1,k+1);
        Uprime1(k+1) = U(1,k+1)-L(1,k+1);
        DprimeN(k+1) = D(N,k+1)+2*U(N,k+1);
        LprimeN(k+1) = L(N,k+1)-U(N,k+1);
    end

    %Boundary
    boundary = [-L(1)*V(1);zeros(N-1,1);-U(N)*V(N+1)];

    %Solve Black-Scholes PDE
    for j=1:M       
        for K=1:N
            %Diagonals for the AI matrix
            AI = zeros(N,N);
            D_AI = [Dprime1(k),D(j,2:N-1),DprimeN(k)];
            U_AI = [Uprime1(k),U(j,2:N)];
            L_AI = [L(j,2:N-1),LprimeN(k)];
            
            AI = diag(L_AI,-1);
            AI = diag(D_AI,0);
            AI = diag(U_AI,1); 
        end
        

    end

    plot(S,V)
    %Use 1D Interpolation to price the option at S0
    %price = interp1(S,V,S0)


    %Tridiagonal Matrix solver
    %Given as part of the course by Dr. Francesco Cosentino
    %Tri-diagonal Matrix Solver

function x = tridiag(D,U,L,B)
  n = length(B);
  x = zeros(n,1); %Vector to store solution

  %Make A upper diagonal
  for j=2:n
      D(j) = D(j) - L(j)*U(j-1) / D(j-1);
      B(j) = B(j) - L(j)*B(j-1) / D(j-1);
  end

  %Back substitution
  x(n) = B(n) / D(n);
  for j=n-1:-1:1
      x(j) = (B(j) - U(j)*x(j+1)) / D(j);
  end
end
