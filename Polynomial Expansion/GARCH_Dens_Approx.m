function [f, M, NM, HM] = GARCH_Dens_Approx(kappa,theta,sigma,X0,N,T)
% Density approximation for a GARCH model
% dX_t = kappa*(theta-X_t)d t + sigma*X_t dW_t

% Input:
% kappa,theta,sigma,X0,T: model parameters
% N: number of moments calculated (2 ≤ N ≤ 4)

% Output:
% M: vector of moments
% NM: vector of normalized moments
% HM: vector of Hermite moments
% f: density approximation up to order N

if(N>4 || N <2)
    error('N is not in [2,3,4]')
end

% Matrix G_N
d1 = (0:N) .* (-kappa + [0 0:N-1] * sigma^2/2);
d2 = (1:N) * kappa * theta;
G = diag(d1, 0) + diag(d2, 1);

% Computing moments
H = X0.^(0:N);
M = H * expm(T*G);
M = M(2:N+1);

% Mean and Variance of X_T
mu = M(2);
sigma = sqrt(M(3)-M(2)^2);

f=0;
NM=0;
HM=0;

end

