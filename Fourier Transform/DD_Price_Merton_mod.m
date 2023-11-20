function [P] = DD_Price_Merton_mod(S,T,r,lambda,alpha,sigma,beta,gamma,L,a,b)

% Input:
% S = Initial price
% r = risk free rate
% lambda,alpha,sigma,beta,gamma = parameters Merton
% L = truncation bound for the integral
% a,b = barriers (already divided by initial price)

% Output: 
% P = Price of a double digital with Maturity T using the characteristic
% function of the price and Fourier pricing rule in the Merton model

% Characteristic function for X_T
CX_T = @(z) exp(T*( 1i*gamma*z - 0.5*sigma^2*z.^2 + lambda*(exp(1i*z*alpha-0.5*beta^2*z.^2)-1) ));

% Characteristic function for DD payoff
CP_T = @(z) ((a/S).^(-1i*z)-(b/S).^(-1i*z))./(1i*z);

% Integrand
I_T = @(z) CP_T(z).*CX_T(z);

% Pricing formula
P = exp(-r*T)/(2*pi)*integral(I_T,-L, L);

end

