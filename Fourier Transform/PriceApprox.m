function [pi] = PriceApprox(N,X0,V0,kappa,sigma,theta,r,rho,v_min,v_max,mu_w,sigma_w,T,k)
% Computes the approximation of the European call option price in the Jacobi model

l = HermiteMoments(N,V0,X0,kappa,sigma,theta,r,rho,v_min,v_max,T,mu_w,sigma_w);
f = FourierCoefficients(N,mu_w,sigma_w,r,T,k);

pi = l*f';

end