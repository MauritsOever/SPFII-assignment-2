function path = GBMsim(drift, sigma, S0, delta_t, T)
%% define variables
N = ceil(T/delta_t);

%% use discretisation of GBM
% dS0 = drift*S0*delta_t + sigma * S0 * dWt = randn
prices = NaN(1,N);
prices(1) = S0;
for i = 2:N
    rnum = randn; 
    dS0 = drift * prices(i-1) * delta_t + sigma*sqrt(delta_t)*prices(i-1)*rnum;
    prices(i) = prices(i-1) + dS0;
end

path = prices;






end
