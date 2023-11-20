function [S] = CIR_Sim(X,T,kappa,theta,sigma,N_time,N_sim,Im)
% Simulates N_sim paths for CIR process 

% Input:
% X: initial value
% T: final time
% kappa, theta, sigma: model parameters
% N_time: number of periods
% N_sim: number of simulations
% Im: 0 if real part only

% Output:
% S: simulated paths

dt  = T/N_time;
S = zeros(N_sim,N_time+1);

% simulation of paths
for i = 1:N_sim
    S(i,1) = X;
    for j=2:N_time+1
        if Im==1
            S(i,j) = S(i,j-1) + kappa*(theta-S(i,j-1))*dt + sigma*sqrt(dt*S(i,j-1))*randn(1);
        else
            S(i,j) = S(i,j-1) + kappa*(theta-S(i,j-1))*dt + sigma*sqrt(dt*max(S(i,j-1),0))*randn(1);
        end
    end
end

end

