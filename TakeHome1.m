% Calculations
r       = 0.1;
sigma   = 0.1;
N       = 2:2:200;
T       = 0.5;
s       = 1;
K       = 0.9;
b       = 1.3;

N_time = 100;
N_sim  = 1e6;

P_B  = zeros(1,length(N));
P_MC = ones(1,length(N)) * MCpriceBarrierUODM(r,sigma,N_time,N_sim,T,s,K,b);

for i = 1:length(N)
    u = 1+r*T/N(i) + sigma*sqrt(T/N(i));
    d = 1+r*T/N(i) - sigma*sqrt(T/N(i));
    P_B(i) = BinomialpriceBarrierUODM(r,d,u,N(i),T,s,K,b);
end

% Plot
hold on
plot(N,P_MC)
plot(N,P_B)
xlabel('N')
ylabel('Option Valuation')
legend('MC Price', 'Binomial Price')
hold off

% We observe that the call option price obtained in the multi-period 
% binomial model with barrier converges (but with oscillation) to the BS 
% call option with barrier price for N goes to infinity.