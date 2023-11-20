function [v] = BSprice(s,K,r,T,sigma)
% Black-Scholes price

t  = 0;
d1 = 1/(sigma*sqrt(T-t))*(log(s/K)+(r+sigma^2/2)*(T-t));
d2 = 1/(sigma*sqrt(T-t))*(log(s/K)+(r-sigma^2/2)*(T-t));

v = s*normcdf(d1)-K*exp(-r*(T-t))*normcdf(d2);

end

