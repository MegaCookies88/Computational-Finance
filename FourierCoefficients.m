function [f] = FourierCoefficients(N,mu_w,sigma_w,r,T,k)
% Computes the first N fourier coefficients for European call option

% Input:
% T: maturity
% k: log strike
% r: risk-free rate

% Output:
% f: coeffs

f = zeros(1,N+1);
I = zeros(1,N+1);

mu = (k-mu_w)/sigma_w;

I(1) = exp(sigma_w^2/2) * normcdf(sigma_w-mu,0,1);            
f(1) = exp(-r*T+mu_w)*I(1) - exp(-r*T+k)*normcdf(-mu,0,1);

for n=1:N
    H = 2^(-0.5*(n-1)) * hermiteH(n-1,mu/sqrt(2));
    I(n+1) = H*exp(mu*sigma_w)*normpdf(mu,0,1) + sigma_w*I(n);
    f(n+1) = exp(-r*T+mu_w)/sqrt(factorial(n))*sigma_w*I(n);
end

end
