%% Ex2

S       = 100;
lambda  = 0.4;
r       = 0;
sigma   = 0.15;
alpha   = -0.5;
beta    = 0.4;
gamma   = r-0.5*sigma^2-lambda*(exp(alpha+0.5*beta^2)-1);
a       = 90;
b       = 100;

P = zeros(4,8);
T = [1/12,1/6,1/4,1/2,1,2,3,4];
L = [10,25,50,100];

for(i = 1:4)
    for(j = 1:8)
        P(i,j) = real(DD_Price_Merton_mod(S,T(j),r,lambda,alpha,sigma,beta,gamma,L(i),a,b)); 
    end
end

%% Ex3

S       = 100;
V       = 0.0175;
r       = 0;
rho     = -0.6;
kappa   = 1.5;
theta   = 0.04;
sigma   = 0.3;
L       = 50;
alpha   = 0.01;

P = zeros(9,8);
Q = zeros(9,8);
T = [1/12,1/6,1/4,1/2,1,2,3,4];
K = [50,80,90,95,100,105,110,120,150];

for(i = 1:9)
    for(j = 1:8)
        [P(i,j),Q(i,j)] = Call_Price_Heston_mod(S,K(i),T(j),r,kappa,theta,sigma,rho,V,alpha,L);
    end
end

% Comparing P and Q
C = max(P-Q,[],'all') % negligible difference

%% Ex3 analysis

S       = 100;
V       = 0.0175;
r       = 0;
rho     = -0.6;
kappa   = 1.5;
theta   = 0.04;
sigma   = 0.3;

P = zeros(9,8);
Q = zeros(9,8);
T = [1/12,1/6,1/4,1/2,1,2,3,4];
K = [50,80,90,95,100,105,110,120,150];

L_tab       = [10,25,50,100];
alpha_tab   = [0.01,0.5,1,1.5,2,5,10];

for(L = L_tab)
    for(alpha = alpha_tab)
        for(i = 1:9)
            for(j = 1:8)
                [P(i,j),Q(i,j)] = Call_Price_Heston_mod(S,K(i),T(j),r,kappa,theta,sigma,rho,V,alpha,L);
                writematrix(round(P,2),"Outputs/HestonPrice_L"+L+"alpha"+alpha+".xls")
            end
        end
    end
end