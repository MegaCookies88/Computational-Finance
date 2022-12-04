%%% Part e %%%

% Defining parameters
S0 = 1;
X0 = log(S0);
sigma = 0.15;
T = 0.5;
lambda = 0.4;
alpha = -0.5;
beta = 0.4;
gamma = -sigma^2/2-lambda*(exp(alpha+0.5*beta^2)-1);
a_min = 0.7;
a_max = 1.3;
L = 50;
nu = -1;

% Defining variables we want to compute
a_tab = linspace(a_min,a_max,100);
orders = 2:30;

P_F = zeros(1,length(a_tab));

F_time = zeros(1,length(a_tab));
CI_time = zeros(1,length(orders));

error = zeros(1,length(orders));

% Computing the Fourier prices
for i=1:length(a_tab)
    tic
    P_F(i) = Digital_Price(a_tab(i),X0,lambda,sigma,alpha,beta,nu,T,L,gamma);
    F_time(i) = toc;
end

% Computing the Interpolated prices
f = @(a) real(Digital_Price(a,X0,lambda,sigma,alpha,beta,nu,T,L,gamma));
for n=1:length(orders)
    tic
    P_CI = ChebInterpol(f,a_tab,orders(n),a_min,a_max);
    CI_time(n) = toc;
    
    error(n) = max(abs(P_F-P_CI));
end

% Plot of log errors
figure(1)
semilogy(orders, exp(-1*orders), 'k--');
hold on
semilogy(orders, error);
grid on
legend('O( e^{-n} )')
ylabel('Absolute interpolation errors')
xlabel('Interpolation order')
title('log maximal absolute error')
hold off

% Plot of computation times
figure(2)
plot(orders, CI_time)
grid on
legend('Interpolation time')
ylabel('Computation time')
xlabel('Interpolation order')
title('Computation time of interpolation for all strikes')
hold off

% Plot of comparison for each strike for order 29
P_CI = ChebInterpol(f,a_tab,29,a_min,a_max);

figure(3)
plot(1:length(a_tab), P_F, 'o')
hold on
plot(1:length(a_tab), P_CI)
legend('Fourier Prices','Interpolated Prices')
ylabel('Price')
xlabel('Strike')
title('Order 29 Interpolated VS Fourier Prices')
hold off



%%% Part d %%%

%  Computes the European Digital Option price in Merton’s model using Fourier
function P = Digital_Price(a,X0,lambda,sigma,alpha,beta,nu,T,L,gamma)

    % Characteristic function of X_T
    CH_X_T = @(z) exp(T*(1i*z*gamma-0.5*sigma^2*z.^2+...
        lambda*(exp(1i*z*alpha-0.5*beta^2*z.^2)-1)));
    
    % Integrand in pricing formula
    integrand = @(u) real(exp(1i*u*X0).* a.^(nu-1i*u)./(1i*u-nu) .* CH_X_T(u+1i*nu));

    % Pricing
    P = exp(-nu*X0)/(pi) * integral(integrand,0,L);

end