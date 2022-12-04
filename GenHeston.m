function [G] = GenHeston(N,r,kappa,theta,sigma,rho,sigma_w)
% Heston model matrix construction

M = 1/2*(N+1)*(N+2);
G = zeros(M,M);

for m = 0:N
    for n = 0:N-m
        ColInd = Ind(m,n); 
        if m > 0
            G(Ind(m-1,n), ColInd) = m*(kappa*theta+sigma^2*(m-1)/2);
        end
        if n > 0
            G(Ind(m,n-1), ColInd) = n*(sigma*rho*m+r)/sigma_w;
            G(Ind(m+1,n-1), ColInd) = -n/(2*sigma_w);
        end
        if n > 1
            G(Ind(m+1,n-2), ColInd) = n*(n-1)/(2*sigma_w^2);
        end
        G(ColInd, ColInd) = -kappa*m;
    end
end

end

