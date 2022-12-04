%% Part B

S0 = 100;
sigma = 0.05;
T = 0.5;
K = 90;
r = 0.01;

N_tab = round(10.^(1:1/3:6));
error = [];

P_BS = blsprice(S0, K, r, T, sigma);

for N=N_tab
    P_MC = BasketMC(S0,sigma,T,K,r,N);
    error = [error abs(P_BS-P_MC)];
end

% plot results
loglog(N_tab, 1.5./sqrt(N_tab), 'k--')
hold on
loglog(N_tab, error)
legend('O(1/sqrt(N))')
ylabel('Absolute MC errors')
xlabel('Number Simulations')
grid on


%% Part D

d = 5;
r = 0.01;
S0 = 100*ones(d,1);
sigma = 0.05*ones(d,1);
T = 0.5;
K_min = 90;
K_max = 110;
N = 1e4;

orders = 2:30;
K_tab = linspace(K_min,K_max,100);

error = zeros(1,length(orders));

P_MC = zeros(1,length(K_tab));
MC_time = zeros(1,length(K_tab));

CI_time = zeros(1,length(orders));

for i=1:length(K_tab)
    tic
    P_MC(i) = BasketMC(S0, sigma, T, K_tab(i), r, N);
    MC_time(i) = toc;
end

for n=1:length(orders)
    tic
    P_CI = ChebInterpol(@(K) BasketMC(S0,sigma,T,K,r,N), K_tab, orders(n), K_min, K_max);
    CI_time(n) = toc;
    
    error(n) = max(abs(P_MC-P_CI));
end

% plot log errors
figure(1)
semilogy(orders, exp(-0.4*orders), 'k--');
hold on
grid on
semilogy(orders, abs(error));
legend('O(exp(-0.5 n))')
ylabel('Absolute interp errors')
xlabel('Interpolation order')
hold off

% plot computation times
figure(2)
plot(orders, CI_time)
grid on
legend('Interpolation time')
ylabel('Computation time')
xlabel('Interpolation order')
hold off

% plot comparison for each strike for order 29
P_CI = ChebInterpol(@(K) BasketMC(S0,sigma,T,K,r,N), K_tab, 29, K_min, K_max);

figure(3)
plot(1:length(K_tab), P_MC, 'o')
hold on
plot(1:length(K_tab), P_CI)
legend('MC Prices','Interpolated Prices')
ylabel('Price')
xlabel('Strike')
hold off






