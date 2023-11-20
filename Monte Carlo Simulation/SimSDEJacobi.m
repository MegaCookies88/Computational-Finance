function [X] = SimSDEJacobi(V0,X0,kappa,sigma,theta,r,rho,v_min,v_max,T,N_time,N_sim)
% Simulates N_sim paths for the Jacobi Process

% Input:
% V0,X0,kappa,sigma,theta,r,rho,v_min,v_max,T: model parameters
% N_time: number of periods
% N_sim: number of simulations

% Output:
% S: simulated paths

dt  = T/N_time;
V = zeros(1,N_sim);
X = zeros(1,N_sim);

% Simulation of paths
for i = 1:N_sim
    V(i) = V0;
    X(i) = X0;
    for j=2:N_time+1
        Z1 = randn(1);
        Z2 = randn(1);
        Q = (V(i)-v_min)*(v_max-V(i))/(sqrt(v_max)-sqrt(v_min))^2;
        X(i) = X(i) + (r-V(i)/2)*dt + rho*sqrt(dt*max(Q,0))*Z1 + sqrt(dt*max(V(i)-rho^2*Q,0))*Z2;
        V(i) = V(i) + kappa*(theta-V(i))*dt + sigma*sqrt(dt*max(Q,0))*Z1;
    end
end

end

