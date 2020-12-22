%% SPF II assignment II, lets get it
% Q1: price atm 1 year spread option: max(ST1-St2-k,0)
% stock 1 is 24-08, so  120$ w/ 17% vol
St1     = 240;
St1_vol = 0.08;

% stock 2 is 17-12, so 80$ w/ 24% vol
St2     = 170;
St2_vol = 0.12;

drift = 0.02; % same for both stocks

% the biggest stock price becomes st1, so the spread stays pos, and k will
% be positive

% correlated with a coefficient of 0.3:
rho = 0.3;
% thinking rn, sim average price of one, then the other, then get payoff,
% then discount...

% the spread option is ATM
% so K = S1 - S2
K = St1 - St2;

% define some params for the sims
N     = 10000; % amount of sims
steps = 250; % assume its 1 year w/ 250 trading days

%% A) block of code that can run some sims/get gap prices

priceA = SpreadOptionMC(St1, St2, St1_vol, St2_vol, K, drift, rho, N, steps);

% seems to work just fine, now to investigate the relationship w/ rho

%% B) investigate the rho relationship
% simming params are defined in code block 1, so run that first

rhos = -1:0.05:1; %vector of different rhos we're gonna try out.
prices_rhos = NaN(1,length(rhos)); %store rho dependend prices

for i = 1:length(rhos)
    prices_rhos(i) = SpreadOptionMC(St1, St2, St1_vol, St2_vol,K,drift,rhos(i),N,steps);
end

plot(rhos, prices_rhos)
% rho up, price down
% makes sense bc theres less extreme differences when rho is high
% so less extremely high payoffs

%% C)
% theoretical question, see Remarkable :DDDD

%% D) spread option price trough GBM:
% not a good idea due to the fact that spread does not follow GBM, but we
% can still calc the price:

% dSpread_t --> drift*spread*dt + sigma_spread * spread * dWt
S0_spread = St1 - St2; % pretend that payoff structure is like a normal option
                          % but S0 = St1 - St2
K_spread = S0_spread; % ATM, so we can compare prices to A)

sigma_spread = sqrt((1*St1_vol)^2 + (-1*St2_vol)^2 + 2* rho * 1 * -1 * St1_vol * St2_vol);

priceD = BlackScholesCall(S0_spread, K_spread, sigma_spread, 1, drift);

% BS is waayyyy lower, because BS assumes the spread cannot become negative
% this is way unrealistic though, spread most definitely can become
% negative

%% E) Bachelier model pricing...Normal version of BS, modified for arithmetic BM





