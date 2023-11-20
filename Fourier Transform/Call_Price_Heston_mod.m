function [P,Q] = Call_Price_Heston_mod(S,K,T,r,kappa,theta,sigma,rho,V,alpha,L)

% Input:
% S = Initial price
% r = risk free rate
% kappa,theta,sigma,rho = parameters Heston
% V = initial vol in Heston model
% L = truncation bound for the integral
% K = strike price
% alpha = damping factor (alpha >0)

% Output:
% P = Price of a European Call with Maturity T and strike K using the characteristic
% function of the price and the Carr-Madan formula - No FFT used
% Q = Price with alternative method subtracting Black Scholes prices

% Characteristic function for Heston model
b = @(nu)(kappa-1i*rho*sigma.*nu);
gamma = @(nu)(sqrt(sigma^2*(nu.^2+1i.*nu)+b(nu).^2));
a = @(nu)(b(nu)./gamma(nu)).*sinh(T*0.5.*gamma(nu));
c = @(nu)(gamma(nu).*coth(0.5*T.*gamma(nu))+b(nu));
d = @(nu)(kappa*theta*T.*b(nu)/sigma^2);
f = @(nu)(1i*(log(S)+r*T).*nu+d(nu));
g = @(nu)(cosh(T* 0.5.*gamma(nu))+a(nu)).^(2*kappa*theta/sigma^2);
h = @(nu)(-(nu.^2+1i.*nu)*V./c(nu));

CHeston_T = @(nu)(exp(f(nu)).*exp(h(nu))./g(nu));

% Characteristic function for BS model
sigmaBS = sqrt(theta);

CBS_T = @(nu)(exp((1i*(log(S)+(r-0.5*sigmaBS^2)*T)*nu)-0.5*(sigmaBS^2)*T*(nu.^2)));

% Integrands
IHeston_T = @(nu)(real((CHeston_T(nu-1i*(alpha+1))./(alpha^2+alpha-nu.^2+1i*(2*alpha+1)*nu)).*exp(-1i*log(K).*nu)));

IBS_T = @(nu) (real((CBS_T(nu-1i*(alpha+1))./(alpha^2+alpha-nu.^2+1i*(2*alpha+1)*nu)).*exp(-1i*log(K).*nu)));
IDiff_T = @(nu)(IHeston_T(nu)-IBS_T(nu));

% Pricing formula
P = (exp(-r*T-alpha*log(K))/pi)*integral(IHeston_T,0,L);

Diff = (exp(-r*T-alpha*log(K))/pi)*integral(IDiff_T,0,L);
Call = BSprice(S, K, r, T, sigmaBS);
Q = Call+Diff;

end

