function [MCprice,MCSE] = monte_carlo(S0,K,B,r,q,T,Nmc,M)
    %Time vector
    Tvals = linspace(0,T,M);

    %Empty matrix to store stock values
    S = zeros(Nmc,M); %Stock prices
    SA = zeros(Nmc,M); %Stock price for antithetic sample
    
    %Empty Matrix to store option payoffs
    hVals = zeros(Nmc,1); 
    hValsAntiThetic = zeros(Nmc,1);
    
    %Generate sample paths using monte-carlo  method
    for i=1:Nmc
        S(i,1) = S0; %Set the initial stock price
        SA(i,1) = S0;
        %Loop through time steps
        for t=1:M-1 
            %Calculate terms here for performance
            Z = randn; %Generate Z from N(0,1)
            sig = vol(S(i,t),Tvals(t));
            DETERM = (r-q-0.5*sig^2)*(Tvals(t+1)-Tvals(t)); %First term in the exponential
            DRIFT = sig*sqrt(Tvals(t+1)-Tvals(t))*Z; %Second term in the exponential

            %Standard sample
            S(i,t+1) = S(i,t)*exp(DETERM+DRIFT);

            %Antithetic sample
            SA(i,t+1) = SA(i,t)*exp(DETERM-DRIFT);

        end
        hVals(i) = payoff(S(i,:),K,B,M);
        hValsAntiThetic(i) = payoff(SA(i,:),K,B,M);
    end
    
    %Compute estimate
    MCprice = (sum(hVals)+sum(hValsAntiThetic))/(2*Nmc);
    


    %Standard error
    subSum = 0;
    NC = Nmc*MCprice; %Calculate this term of the sum here for efficiency
    for i = 1:Nmc
        subSum = subSum + ((hVals(i)+hValsAntiThetic(i))/2)^2;
    end
    MCSE = sqrt((subSum-NC)/(Nmc*(Nmc-1)));
end

%Payoff for up-and-out barrier option
function h = payoff(Svals,K,B,M)
    if max(Svals) < B
        h = max(Svals(M)-K,0);
    else
        h = 0;
    end
end

%Function used here for calculating volatility
function sig = vol(s,t)
    alpha = 0.35;
    sig = 0.25*exp(-t)*(100/s)^alpha;
end