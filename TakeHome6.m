%% PREPARING DATA

training_size = 2500;
testing_size = 5000;

%%% Preparing a training set %%%
% Parameters intervals and sampling for training
strike = 0.4 + (1.6-0.4).*rand(training_size,1);
T = 11/12 + (1-11/12).*rand(training_size,1);
r = 0.015 + (0.025-0.015).*rand(training_size,1);
kappa = 1.4 + (2.6-1.4).*rand(training_size,1);
theta = 0.45 + (0.75-0.45).*rand(training_size,1);
rho = -0.75 + (-0.45+0.75).*rand(training_size,1);
sigma = 0.01 + (0.1-0.01).*rand(training_size,1);
v0 = 0.01 + (0.1-0.01).*rand(training_size,1);

% Training Set
Heston_train = zeros(training_size,8);
Heston_train(:,1) = strike(randperm(training_size));
Heston_train(:,2) = T(randperm(training_size));
Heston_train(:,3) = r(randperm(training_size));
Heston_train(:,4) = kappa(randperm(training_size));
Heston_train(:,5) = theta(randperm(training_size));
Heston_train(:,6) = rho(randperm(training_size));
Heston_train(:,7) = sigma(randperm(training_size));
Heston_train(:,8) = v0(randperm(training_size));

% Computing training prices using FFT
Heston_price_train = zeros(training_size,1);
for i = 1:training_size
    strike = Heston_train(i,1);
    T = Heston_train(i,2);
    r = Heston_train(i,3);
    kappa = Heston_train(i,4);
    theta = Heston_train(i,5);
    rho = Heston_train(i,6);
    sigma = Heston_train(i,7);
    v0 = Heston_train(i,8);
    S = 1;
    alpha = 2;
    Heston_price_train(i,1) = FFT_Heston(kappa,theta,sigma,rho,r,v0,S,strike,T,alpha);
end

%%% Preparing a testing set %%%
% Parameters intervals and sampling for testing
strike = 0.5 + (1.5-0.5).*rand(testing_size,1);
T = 11/12 + (1-11/12).*rand(testing_size,1);
r = 0.015 + (0.025-0.015).*rand(testing_size,1);
kappa = 1.5 + (2.5-1.5).*rand(testing_size,1);
theta = 0.5 + (0.7-0.5).*rand(testing_size,1);
rho = -0.7 + (-0.5+0.7).*rand(testing_size,1);
sigma = 0.02 + (0.1-0.02).*rand(testing_size,1);
v0 = 0.02 + (0.1-0.02).*rand(testing_size,1);

% Testing Set
Heston_test = zeros(testing_size,8);
Heston_test(:,1) = strike(randperm(testing_size));
Heston_test(:,2) = T(randperm(testing_size));
Heston_test(:,3) = r(randperm(testing_size));
Heston_test(:,4) = kappa(randperm(testing_size));
Heston_test(:,5) = theta(randperm(testing_size));
Heston_test(:,6) = rho(randperm(testing_size));
Heston_test(:,7) = sigma(randperm(testing_size));
Heston_test(:,8) = v0(randperm(testing_size));

%%% Fitting a GPR model from MATLAB then predicting %%%
gprMdl1 = fitrgp(Heston_train,Heston_price_train,'KernelFunction','squaredexponential');
gprMdl1_predict_Heston = predict(gprMdl1,Heston_test);

%% FITGPR homemade

bounds = [1,2;0.3,0.6;0.001,0.01];

% squared exponential kernel
theta0 = [2,0.47,0.0014];
[ML,HP,PM,PC] = fitGPR(Heston_train,Heston_price_train,@sqrdexp,theta0,bounds,Heston_test);
Optimal_Parameters = HP

%% PLOT

plot(gprMdl1_predict_Heston,PM)
hold on
%plot([1,training_size],[1,training_size])
xlabel('Heston Matlab GPR')
ylabel('Heston Homemade GPR')
title('Comparing to Matlab GPR')

mean_error = mean(abs(gprMdl1_predict_Heston-PM))

%% Kernels

%%% Kernels functions and gradients %%%
% Squared Exponential
function [cov,grad_s,grad_l] = sqrdexp(Xn,Xm,sigma,length_scale)
    l = length_scale;
    z = (pdist2(Xn,Xm).^2)/(2*l^2);
    cov = (sigma^2)*exp(-z);
    grad_s = 2*sigma*exp(-z);
    grad_l = (sigma^2).*(pdist2(Xn,Xm).^2)/(l^3).*exp(-z);
end

% Linear
function [cov,grad_s0,grad_s1,grad_c] = linear_kernel(Xn,Xm,sigma0,sigma1,c)
    cov = (sigma0^2)+(sigma1^2)*(Xn-c)*(Xm-c)';
    grad_s0 = 2*sigma0;
    grad_s1 = 2*sigma1;
    grad_c  = (sigma1^2)-(Xm+Xn)+2*c;
end

% Periodic
function [cov,grad_s,grad_l,grad_p] = periodic_kernel(Xn,Xm,sigma,length_scale,p)
    l = length_scale;
    z = pi*pdist2(Xn,Xm)/p;
    cov = (sigma^2)*exp(-(2/(l^2))*(sin(z).^2));
    grad_s = 2*sigma*exp(-(2/(l^2))*(sin(z).^2));
    grad_l = (sigma^2)*((4/(l^3))*(sin(z).^2))*exp(-(2/(l^2))*(sin(z).^2));
    grad_p = (sigma^2)*((4/(l^2))*(z/p)*sin(z)*cos(z))*exp(-(2/(l^2))*(sin(z).^2));
end