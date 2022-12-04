% Part F

N = 50;          
X0 = 0;
V0 = 0.04;
kappa = 0.5;
theta = 0.04;
sigma = 1;
r = 0;
rho = -0.5;
T = 1/12;
v_min = 1e-4;
v_max = 0.08;

mu_w = X0 + (r-1/2*sigma^2)*T;                     
sigma_w = sigma^2*T;

N_sim  = 1e6;
N_time = 100;

k = -0.1;

% polynomial expansion approach
tic
PE_value = PriceApprox(N,X0,V0,kappa,sigma,theta,r,rho,v_min,v_max,mu_w,sigma_w,T,k);
PE_t = toc;

% standard Monte Carlo approach
tic
X = SimSDEJacobi(V0,X0,kappa,sigma,theta,r,rho,v_min,v_max,T,N_time,N_sim);
MC_value = exp(-r*T) * mean(exp(-r*T)*max(exp(X)-exp(k),0));
MC_t = toc;

% comparing absolute error
err = abs(PE_value - MC_value);
formatSpec = 'The absolute error is %1.2d \n';
fprintf(formatSpec,err)

% comparing computation time
formatSpec = 'PE takes %4.4f seconds and MC takes %4.4f seconds \n';
fprintf(formatSpec,PE_t,MC_t)

