function [ML,HP,PM,PC] = fitGPR(X,y,kernel,theta0,bounds,X_test)
% Fits a Gaussian process regression

% Inputs:
% X/X_test: training/testing input
% y/y_test: training/testing output
% kernel: kernel used
% params: kernel parameters
% bounds: parameters bounds

% Outputs:
% ML: maxima of marginal likelihood
% HP: optimal hyper-parameters
% PM: posterior mean
% PC: posterior covariance

% optimization function 
f = @(x) -log_marginal_likelihood(X,y,kernel,x);
% initial values for optimization
x0 = theta0; 
% boundaries for variables
LB = bounds(:,1);
UB = bounds(:,2);
% likelihood maximization
HP = fminsearchcon(f,x0,LB,UB);
ML = log_marginal_likelihood(X,y,kernel,HP);

% use optimized model
N = size(X,1);
[K] = kernel(X,X,HP(1),HP(2));
[K_star] = kernel(X_test,X,HP(1),HP(2));
[K_star2] = kernel(X_test,X_test,HP(1),HP(2));
K_y = K + (HP(3)^2)*eye(N);

% I don't use this because it doesn't work it is not stable and gives huge values
%L = chol(K_y);
%alpha = L'\(L\y);
%v = L\K_star;
%PM = K_star*alpha;
%PC = K_star2 - v'*v ;

PM = K_star*(K_y\y);
PC = K_star2-K_star*(K_y\(K_star'));

end

% log marginal likelihood
function [l] = log_marginal_likelihood(X,y,kernel,theta)
    N = size(X,1);
    [K] = kernel(X,X,theta(1),theta(2));
    K_y = K + (theta(3)^2)*eye(N);
    alpha = K_y\y;
    
    l = (-1/2)*(y'*alpha + log(det(K_y))/2 + N*log(2*pi)/2);
end