function BS = BlackScholesCall(S0, K, sigma, T, r)

d1 = (log(S0/K) + (r+(sigma^2)/2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);

BS = S0*normcdf(d1) - K*exp(-r*T)*normcdf(d2);

end