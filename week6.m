%% Ex1

X = 1;
kappa  = 0.5;
theta  = 0.4;
sigma  = 0.2;
T = 0.5;
N = 4;

[f, M, NM, HM] = GARCH_Dens_Approx(kappa,theta,sigma,X,N,T);
round(M,2)

%% Ex2

X0 = 5.1;
V0 = 0.04;
kappa  = 1;
theta  = 0.04;
sigma  = 0.2;
r = 0.03;
rho = -0.8;
T = 1/52;

N = 4;

% Assume mu_w=0 and sigma_w=1 to compute mu_w and sigma_w ...
sigma_w = 1;

G = GenHeston(2,r,kappa,theta,sigma,rho,sigma_w);

H = [1, V0, X0, V0^2, V0*X0, X0^2];
moments = H * expm(G*T);

mu_w = moments(3);
sigma_w = sqrt(moments(6) - mu_w^2);

% We now compute the 4 moments as asked in the exercise
G = GenHeston(N,r,kappa,theta,sigma,rho,sigma_w);

X0_tilde = (X0-mu_w) / sigma_w;
H = zeros(length(G),1);

for m = 0:N
    for n = 0:N-m
        ColInd = Ind(m,n);
        H(ColInd, 1) = V0^m * X0_tilde^n;
    end
end

moments = H' * expm(G*T);

% We print the needed moments
fprintf('First moment: %.4f\n', moments(Ind(0,1)))
fprintf('Second moment: %.4f\n', moments(Ind(0,2)))
fprintf('Third moment: %.4f\n', moments(Ind(0,3)))
fprintf('Fourth moment: %.4f\n', moments(Ind(0,4)))

%% Draft

% Generalized Hermite polynomials
NM = zeros(N+1,1);
for n=1:N+1
    for m=0:n-1
        NM(n) = NM(n) + (nchoosek(n-1,m)*(-mu_w)^(n-1-m)/(sigma_w)^(n-1))*M(m+1);
    end
end

syms y
HM = zeros(N+1,1);
for n=1:N+1
    for m=0:n-1
        hermite_coeff = coeffs(hermiteH(n-1,y),'ALL');
        HM(n) = 2^(-0.5*(n-1))/sqrt(factorial(n-1))*(hermite_coeff * NM(n:-1:1));
    end
end