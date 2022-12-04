training_size = 1000;
testing_size = 500;

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
tic
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
train_FFT_time = toc;


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

% Computing testing prices using FFT
Heston_price_test = zeros(testing_size,1);
tic
for i = 1:testing_size
    strike = Heston_test(i,1);
    T = Heston_test(i,2);
    r = Heston_test(i,3);
    kappa = Heston_test(i,4);
    theta = Heston_test(i,5);
    rho = Heston_test(i,6);
    sigma = Heston_test(i,7);
    v0 = Heston_test(i,8);
    S = 1;
    alpha = 2;
    Heston_price_test(i,1) = FFT_Heston(kappa,theta,sigma,rho,r,v0,S,strike,T,alpha);
end
test_FFT_time = toc;


%%% Fitting a GPR model %%%
tic
gprMdl1 = fitrgp(Heston_train,Heston_price_train,'KernelFunction','squaredexponential');
GPR_time = toc;


%%% Predicting and performance analysis
tic
predict_Heston=predict(gprMdl1,Heston_test);
pred_time = toc;

plot(predict_Heston,Heston_price_test)
title('Out of Sample Prediction')
xlabel('Heston GPR')
ylabel('Heston FFT')
mae_error = max(abs(predict_Heston-Heston_price_test));
aae_error = mean(abs(predict_Heston-Heston_price_test));
