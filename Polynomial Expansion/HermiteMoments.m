function [l] = HermiteMoments(N,V0,X0,kappa,sigma,theta,r,rho,v_min,v_max,T,mu_w,sigma_w)
% Computes the first N Hermite moments

% Input:
% V0,X0,kappa,sigma,theta,r,rho,v_min,v_max,T,mu_w,sigma_w: model parameters
% N: number of moments calculated

% Output:
% l: Hermite moments

% Position in basis vector
Ind = @(m,n) (m+n)*(m+n+1)/2 + 1 + n;

% Hermite polynomials
H = @(x,n) 2^(-0.5*n) * hermiteH(n,x/sqrt(2));

% Basis vector B_N(V_0,X_0)
B = zeros(1,(N+2)*(N+1)/2);

for m = 0:N
    for n = 0:N-m
        H_n = H((X0-mu_w)/sigma_w,n)/sqrt(factorial(n)); % H_n(X_0)
        B(1,Ind(m,n)) = V0^m * H_n;
    end
end

% Matrix G_N
M = 1/2*(N+1)*(N+2);
G = zeros(M,M);

for m = 0:N
    for n = 0:N-m
        ColInd = Ind(m,n); 
        G(ColInd,ColInd) = -kappa*m-(sigma^2*m*(m-1))/(2*(sqrt(v_max)-sqrt(v_min))^2);
        if(m>=2)
            G(Ind(m-2,n),ColInd)=...
            -(sigma^2*m*(m-1)*v_max*v_min)/(2*(sqrt(v_max)-sqrt(v_min))^2);
        end
        if(m>=1 && n>=1)
            G(Ind(m-1,n-1),ColInd)=...
            -(sigma*rho*m*sqrt(n)*v_max*v_min)/(sigma_w*(sqrt(v_max)-sqrt(v_min))^2);
        end
        if(m>=1)
            G(Ind(m-1,n),ColInd)=...
            kappa*theta*m+(sigma^2*m*(m-1)*(v_max+v_min))/(2*(sqrt(v_max)-sqrt(v_min))^2);
        end
        if(n>=1)
            G(Ind(m,n-1),ColInd)=...
            r*sqrt(n)/sigma_w + (sigma*rho*m*sqrt(n)*(v_max+v_min))/(sigma_w*(sqrt(v_max)-sqrt(v_min))^2);
        end
        if(n>=2)
            G(Ind(m+1,n-2),ColInd)=...
            sqrt((n)*(n-1))/(2*sigma_w^2);
        end
        if(n>=1)
            G(Ind(m+1,n-1),ColInd)=...
            -sqrt(n)/(2*sigma_w)-(sigma*rho*m*sqrt(n))/(sigma_w*(sqrt(v_max)-sqrt(v_min))^2); 
        end 
    end
end

%l = B * expm(G*T);
%l = l(1:N+1);

M = eye((N+2)*(N+1)/2);
l = zeros(1,N+1);
for n=0:N
    l(n+1) = B * expm(G*T) * M(:,Ind(0,n));
end

end