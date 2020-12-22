function price = SpreadOptionMC(St1, St2, St1_vol, St2_vol, K, drift, rho, N, steps)
%% function to get spread option prices for two stocks that follow GBM
% their drift is assumed to be the same under Q (risk-free rate)

%% variable explanation
% St1 = price of first stock
% St2 = prices of second stock (St1>St2)
% St1_vol is (yearly) volatility of stock 1
% St2_vol is (yearly) volatility of stock 2
% K is strike of the spread option (ATM would be K = St1-St2
% drift is assumed risk free rate that the stocks make in expectation
% rho is the correlation between the stocks' driving brownian motions, as
% they are assumed to follow seperate GBM's

% N is amount of simulated price paths
% steps is amount of steps in a year (as we assume a fixed maturity of a year)

%% actual code part

endpricesS1 = NaN(1,N);
endpricesS2 = NaN(1,N); % vector to store endprices in

payoff = NaN(1,N); % vector to store payoffs in

for i = 1:N % for every sim
    pricesS1 = NaN(1,steps); % clear/create vectors that store prices
    pricesS2 = NaN(1,steps); % same
    
    pricesS1(1) = St1;
    pricesS2(1) = St2; % initialize
    
    for t = 2:steps % for every step
        epS1 = randn;
        epS2 = rho*epS1 + randn*sqrt(1-rho^2); % create correlated random normal samples
        
        % calculate differential step, assume GBM
        % dSt = mu*st*dt + sigma*sqrt(dt)*st*eps
        dS1 = drift*pricesS1(t-1)*(1/steps) + St1_vol*sqrt(1/steps)*pricesS1(t-1)*epS1;
        dS2 = drift*pricesS2(t-1)*(1/steps) + St2_vol*sqrt(1/steps)*pricesS2(t-1)*epS2;
        
        % add to price(t-1) so you get price(t)
        pricesS1(t) = pricesS1(t-1) + dS1;
        pricesS2(t) = pricesS2(t-1) + dS2;
        
        % "back, jack, and do it again" - steely dan
        % such a great song
    end
    
    % store endprices in the endprices vector for each simmed path
    endpricesS1(i) = pricesS1(length(pricesS1)); % apparently faster then pricesS1(end)
    endpricesS2(i) = pricesS2(length(pricesS2)); 
    
    % calculate payoff and store in payoff vector for each sim
    payoff(i) = max(endpricesS1(i)-endpricesS2(i)-K,0); % payoff structure for spread option
    
end


% calc option price as average i.e. expected payoff of simulated paths:
price = mean(payoff);




end