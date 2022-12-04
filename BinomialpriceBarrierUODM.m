function [P] = BinomialpriceBarrierUODM(r,d,u,N,T,s,K,b)
% Input:
% r   = risk free interest rate
% d,u = expected returns on the underlying asset
% N   = number of time-intervals
% T   = maturity
% s   = initial asset price
% K   = strike price
% b   = barrier - threshold 
% Output:
% P = price of an European option in the multi-period binomial model 

if mod(N,2) ~= 0
    error('N must be even.')
end

if  (s < 0 || b <= 0 || K < 0)   
    error('Prices must be positive.')
end

rt  = r * T/N;
pu  = (1+rt-d)/(u-d);
pq  = 1-pu;

C = zeros(N+1,N+1);

% valuation of payoff at maturity nodes given the barrier b
for i = 0:N/2
    for j = 0:N/2
        S_half = s * (u^(i)) * (d^(N/2-i)); % price at T/2 
        S = S_half * (u^(j)) * (d^(N/2-j)); % price at T
        if (S < b && S_half < b)
            C(i+j+1) = max(S - K, 0);
        end
    end
end

% valuation of payoff at previous nodes
for k = N-1:-1:0
    for i = 0:k
        C(i+1) = 1/(1+rt) * (pu*C(i+2)+pq*C(i+1));
    end
end

% computing binomial price
P = C(1);

end

