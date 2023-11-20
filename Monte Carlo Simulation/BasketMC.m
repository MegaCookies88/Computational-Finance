function [P] = BasketMC(S0,sigma,T,K,r,N)
% Computes basket option price using MC 

% Input:
% S0,sigma,T,K,r: model parameters
% N: number of simulations

% Output:
% P: Option price

d = length(S0);
sigma_hat = repmat(sigma',N,1);

W = sqrt(T).*sigma_hat.*randn(N,d);
paths = repmat(S0',N,1).*exp(T*(r-0.5*sigma_hat.^2)+W);
paths = sum(paths,2)./d;

payoffs = max(max(paths,[],2)-K,0);

P = exp(-r*T) * mean(payoffs);

end

