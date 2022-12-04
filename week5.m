%% Ex1

X      = 0.04;
T      = 1/12;
kappa  = 0.5;
theta  = 0.04;
sigma  = 2;
N_sim  = 100;
N_time = 100;

CIR_S = CIR_Sim(X,T,kappa,theta,sigma,N_time,N_sim,1);

CIR_S_Re = CIR_Sim(X,T,kappa,theta,sigma,N_time,N_sim,0);

plot(CIR_S_Re(1:10,:)')
legend("Sim "+string(1:10))

%% Ex2

M = zeros(1,3);

G = GenCIR(kappa,theta,sigma,3);

for n=1:3
    p = zeros(1,4);
    p(n+1) = 1;
    M(n) = round(MonCIR(G,X,T,p'),4);
end

N_sim  = 1e6;
N_time = 100;

M_S = zeros(1,3);
CIR_S_Re = CIR_Sim(X,T,kappa,theta,sigma,N_time,N_sim,0);

for n=1:3
    M_S(n) = round(mean(CIR_S_Re(:,N_time+1).^n),4);
end
