function [G] = GenCIR(kappa,theta,sigma,N)
% Computes generator matrix for CIR process 

% Input:
% kappa, theta, sigma: model parameters
% N: polynomial max dimension 

% Output:
% G: generator matrix

d1 = (0:N)* (-kappa);
d2 = (1:N)* kappa * theta + 0.5 * sigma^2 * (1:N).* (0:N-1);

G = diag(d1, 0) + diag(d2, 1);

end

