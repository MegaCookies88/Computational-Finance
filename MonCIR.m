function [M] = MonCIR(G,X,T,p)
% Computes moments for CIR process 

% Input:
% G: generator matrix
% X: initial value
% T: final time
% p: polynomial coordinate

% Output:
% M: moment

N = size(G);
N = N(1)-1;

H_N  = X .^ (0:N); % canonical basis of polynomials at X

%M = H_N * expmv(T,G,p);
M = H_N * expm(G*T) * p;

end

