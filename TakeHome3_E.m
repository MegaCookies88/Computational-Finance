% Part E

N = 10;          
X0 = 0;
V0 = 0.04;
kappa = 0.5;
theta = 0.04;
sigma = 1;
r = 0;
rho = -0.5;
T = 1/12;
v_min = 1e-4;
v_max = 0.08;

mu_w = X0 + (r-1/2*sigma^2)*T;                     
sigma_w = sigma^2*T;

L = 100;
alpha = 1;

log_strike = linspace(-0.1,0.1,50);
IV = zeros(2,50);

for k=1:50
    value_Hermite = PriceApprox(N,X0,V0,kappa,sigma,theta,r,rho,v_min,v_max,mu_w,sigma_w,T,log_strike(k));
    value_Heston  = Call_Price_Heston_mod(exp(X0),exp(log_strike(k)),T,r,kappa,theta,sigma,rho,V0,alpha,L);
    IV(1,k) = blsimpv(exp(X0), exp(log_strike(k)), r, T, value_Hermite);
    IV(2,k) = blsimpv(exp(X0), exp(log_strike(k)), r, T, value_Heston);
end

plot(log_strike,100*IV)
ylabel('Implied Volatility in %')
xlabel('Log Strike')
title('Implied volatility smile for the Hermite and Heston models')
legend('Hermite model','Heston model','Location','northeast')