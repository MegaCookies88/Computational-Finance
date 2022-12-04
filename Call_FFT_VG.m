function [P,beta,tmp] = Call_FFT_VG(K, T, r, theta, gamma, sigma, S, alpha, rule)
% Computes the value of call option option using FFT technique for a VG model

% Input:
% K: strike price of the call option
% T: maturity of the call option
% r: risk free rate
% theta, gamma, sigma: parameters of the VG model
% S: initial stock price
% alpha: damping factor
% rule: weights rule (0: trapezoid, 1: simpson)

% Output:
% P: price of the call option with strike K and maturity T
% beta: discretization log strikes
% tmp: array of prices produced by algorithm

% Charasteristic function for X_T of VG model
b = @(u)(log(1-1i*theta*gamma.*u+sigma^2*0.5*gamma.*u.^2));
a = @(u)(exp(-(T/gamma)*b(u)));

CX_T = @(u) (exp(1i*u.*(r*T-log(a(-1i)))).*a(u));

% FFT setup
eta = 0.25;
N = 2^12;
lambda = 2*pi / (N * eta);

beta = -(N-1)*lambda/2+lambda.*(0:N-1); %discretization of log-strike 
nu = eta .* (0:N-1); %discretization of integration variable

logStrikes = log(K/S);

% Weights integration
trapezoid_weights = ones(1,N);
trapezoid_weights(1) = 0.5;
trapezoid_weights(N) =0.5;

Simpson_weights = repmat([2,4],1,N/2);
Simpson_weights(1) = 1;
Simpson_weights = Simpson_weights/3; 

if rule==0
    weights = trapezoid_weights;
else
    weights = Simpson_weights;
end

% Linear interpolation
ind = min(max(1+floor( (logStrikes - beta(1)) / lambda ),1),N-1);
diff = logStrikes - beta(ind);
nu_tmp = (nu - 1i*(alpha+1)); % Points where ch. function will be evaluated

cgf2 = CX_T(nu_tmp); % Evaluation of ch. function at nu_tmp

% FFT Algorithm
func = exp(-1i*beta(1)*nu) .* cgf2 ./ ((alpha + 1i*nu).*(1+alpha+1i*nu)) .* eta .* weights;
tmp = (exp(-alpha.*beta-r*T) ./ pi) .* real(fft(func));

% Pricing formula
P = S*(tmp(ind) + diff .* (tmp(ind+1) - tmp(ind)) ./ lambda);

end

