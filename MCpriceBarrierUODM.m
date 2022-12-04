function [P] = MCpriceBarrierUODM(r,sigma,N_time,N_sim,T,s,K,b)
% Input:
% r      = risk free interest rate
% sigma  = standard deviation
% N_time = number of time-intervals
% N_sim  = number of time-intervals
% T      = maturity
% s      = initial asset price
% K      = strike price
% b      = barrier - threshold 
% Output:
% P = MC price

if mod(N_time,2) ~= 0
    error('N_time must be even.')
end

if  (s < 0 || b <= 0 || K < 0)   
    error('Prices must be positive.')
end

dt  = T/N_time;

S = zeros(N_sim,N_time+1);
C = zeros(1,N_sim);

% simulation of paths
for i = 1:N_sim
    S(i,1) = s;
    for j=2:N_time+1
        S(i,j) = (1 + r*dt + sigma*sqrt(dt)*randn(1)) * S(i,j-1);
    end
end

% valuation of payoff
for i = 1:N_sim
    if (S(i,N_time+1) < b && S(i,N_time/2+1) < b)
        C(i) = max(S(i,N_time+1) - K, 0);
    end
end

% computing price MC
P = exp(-r*T) * mean(C);

end