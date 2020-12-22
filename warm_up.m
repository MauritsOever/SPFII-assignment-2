%% 1) warm up
% in accordance w/ the assignment warm up protocol
% got to sim a SABR model, stochastic vol
% EZ enough right

%% define some variables for SABR-tooth model
% F is underlying, defined as forward rate
% dFt = sigmat * (Ft^beta) * dW1(t)
% dSigma = alpha * sigmat * dW2(t)
% wieners are correlated by rho

% lets define some values:
F0 = 1;
r = 0; % unused i guess
sigma_0 = 0.2;
beta = 1;
rho = -0.7;
alpha = 0.6;

steps = 250; % amount of time steps

% generating correlated random samples
% eps1 = randn
% eps2 = rho*eps1 + randn2*sqrt(1-rho^2)

%% code chunk that can get price path of F:
F = NaN(1,steps); % get vector to store F path, preallocation for speed ofc
F(1) = F0; % initialize

sigma = NaN(1,steps); % get vector to store sigmat path
sigma(1) = sigma_0; % initialize

for t = 2:steps % generate a single path with steps amount of steps, here it is daily
    eps1 = randn;
    eps2 = rho*(1/steps) * eps1 + randn * sqrt(1-(rho*(1/steps))^2); % get randns w/ correlation
    
    dFt     = sigma(t-1) * F(t-1)^beta * eps1 * sqrt(1/steps);
    dSigmat = exp(-0.5*alpha^2*(1/steps)+alpha*eps2*sqrt(1/steps)); % calc changes per t
    
    F(t)     = F(t-1) + dFt;
    sigma(t) = sigma(t-1) * dSigmat; % add changes to initials for next observation
    
    % do that (steps) times, assumed this is daily again
end

plot(sigma) % ayy lmao it works fam


%% plot hist of returns and compare to normal dist with some mean and variance
N = 100; % amount of sims

returns = []; % vector that can store returns, get a lot of them 
                              % this method sims 100 years of daily returns
                              % so should be pretty stable on average.

for i = 1:N
    F = NaN(1,steps); % get vector to store F path, preallocation for speed ofc
    F(1) = F0; % initialize

    sigma = NaN(1,steps); % get vector to store sigmat path
    sigma(1) = sigma_0; % initialize

    for t = 2:steps % generate a single path with steps amount of steps, here it is daily
        eps1 = randn;
        eps2 = rho*(1/steps) * eps1 + randn * sqrt(1-(rho*(1/steps))^2); % get randns w/ correlation
    
        dFt     = sigma(t-1) * F(t-1)^beta * eps1 * sqrt(1/steps);
        dSigmat = exp(-0.5*alpha^2*(1/steps)+alpha*eps2*sqrt(1/steps)); % calc changes per t
    
        F(t)     = F(t-1) + dFt;
        sigma(t) = sigma(t-1) * dSigmat; % add changes to initials for next observation
    
        % do that (steps) times, assumed this is daily again
    end
    
    % convert path to log rets and append to storage vector
    returns = [returns,price2ret(F)];
end


histfit(returns, 150 ,'Normal') % fit a normal pdf to the histogram of returns
xlim([-0.1 0.1]) % thick tails, leptokurtic distribution, as expected since we're talking about returns
% end of warm up 
% or is it xdd

%% Get sims of SABR, average payoffs and get discount (option pricing MC) - ATM first
% define some variables for SABR-tooth model:

% F is underlying, defined as forward rate
% dFt = sigmat * (Ft^beta) * dW1(t)
% dSigma = alpha * sigmat * dW2(t)
% wieners are correlated by rho

% lets define some values:
F0 = 1;
r = 0; % unused i guess
sigma_0 = 0.2;
beta = 1;
rho = -0.7;
alpha = 0.6;

T = 1; % one year
steps = 250; % amount of time steps, or T/delta_t

moneyness = 0; % moneyness in percent of S0/F0...
K = F0 - moneyness*F0;

N = 100; % amount of sims

payoff = []; % vector that can store returns, get a lot of them 
                              % this method sims 100 years of daily returns
                              % so should be pretty stable on average.

for i = 1:N
    F = NaN(1,steps); % get vector to store F path, preallocation for speed ofc
    F(1) = F0; % initialize

    sigma = NaN(1,steps); % get vector to store sigmat path
    sigma(1) = sigma_0; % initialize

    for t = 2:steps % generate a single path with steps amount of steps, here it is daily
        eps1 = randn;
        eps2 = rho*(1/steps) * eps1 + randn * sqrt(1-(rho*(1/steps))^2); % get randns w/ correlation
    
        dFt     = sigma(t-1) * F(t-1)^beta * eps1 * sqrt(1/steps);
        dSigmat = exp(-0.5*alpha^2*(1/steps)+alpha*eps2*sqrt(1/steps)); % calc changes per t
    
        F(t)     = F(t-1) + dFt;
        sigma(t) = sigma(t-1) * dSigmat; % add changes to initials for next observation
    
        % do that (steps) times, assumed this is daily again
    end
    
    % convert path to log rets and append to storage vector
    payoff(i) = max(F(length(F)) - K,0);
    mean_payoff = mean(payoff);
    SIMprice = mean_payoff * exp(-r*1); % mat = 1 year
end


% BS = BlackScholesCall(S0, K, sigma, T, r)
% we need to figure out sigma in a way...

close = SIMprice - BlackScholesCall(F0, K, sigma_0, 1, r);
                            % taking intial vol already gets you pretty close, 
                            % but lets try average sigma
close2 = SIMprice - BlackScholesCall(F0, K, mean(sigma), 1, r);
                            % similar results, maybe the numbers are too
                            % small, lets try geometric average
close3 = SIMprice - BlackScholesCall(F0, K, geomean(sigma), 1, r);
                            % initial vol seems to be the best...
                            % makes sense since your pricing at t=0...
close4 = SIMprice - BlackScholesCall(F0, K, std(F), 1, r); 
                            % use realized vol over F0??
                            % worst so far... darnit..the search
                            % continues..

close5 = SIMprice - BlackScholesCall(F0, K, blsimpv(F0, K, r, 1, SIMprice), 1, r);
                            % can get implied vol as vol
                            % yeet, very close, but probably not the vol
                            % we're looking for 

