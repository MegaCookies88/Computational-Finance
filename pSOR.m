function [xnew,k,errest] = pSOR(A,b,omega,phi,tol,itermax,xinit)


%function [x,k,errest] = pSOR(A,b,omega,phi,tol,itermax,xinit)
%
% solves the linear complementarity system 
%
% Ax-b > 0 
% x-phi> 0
% (Ax-b)'(x-phi)=0
%
% by pSOR method (with parameter omega), xinit being the initial guess
%
% outputs:
%
% x : solution of the pSOR algorithm
% k : number of iterations performed 
% errest : error estimate computed as the relative increment 



% The pSOR algorithm should stop as soon as the relative difference between two consecutive 
% iterations (variable errest) is lower than a prescribed tolerance, 
% or if the maximum number of iterations (variable k) per time-step has been reached.
% Xold, errest and k should be properly initialized

errest=1;
k=0;
xold=xinit;
xnew=0*xold;
N=length(xold);

while errest>tol && k<itermax
    k=k+1;
    for i = 1:N
        % note that sum_{j=1}^{i-1} a_ij x_new(j) can be written as a row-column vector multiplication
        t = omega/A(i,i)*( b(i) - A(i,1:i-1)*xnew(1:i-1) - A(i,i+1:end)*xold(i+1:end) ) + (1-omega)*xold(i);
        xnew(i) = max(t,phi(i));
    end
    errest=norm(xnew-xold)/norm(xold);
    xold=xnew;
end
