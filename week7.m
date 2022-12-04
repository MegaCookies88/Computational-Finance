%% Ex1

exact_sol = @(x) 3*x.^2 + 1;

sigma = 1;
r = 1;
x_min = 0;
x_max = 4;
Nx = 251;

rhs = @(x) 3*r*x.^2 -(r-sigma^2/2)*6*x+(r-3*sigma^2); % forcing

bc_left = exact_sol(x_min);
bc_right = exact_sol(x_max);

[V,FD_grid] = bs_logform_steady(sigma,r,rhs,bc_left,bc_right,x_min,x_max,Nx);

% comparing solutions
figure
set(gca,'FontSize',20)
plot(FD_grid,exact_sol(FD_grid),'ob','LineWidth',2,'DisplayName','Exact solution') 
hold on
plot(FD_grid,V,'r','LineWidth',2,'DisplayName','Finite diff. approx.')
 
grid on
legend show
set(legend,'Location','NorthWest')

% difference between exact solution and finite differences approximation
max(abs(V'-exact_sol(FD_grid)))

%% Ex2

exact_sol = @(x,t) 4*t*(3*x.^2 + 1);

sigma = 1;
r = 1;
x_min = 0;
x_max = 1;
T = 1;
Nx = 21;

forcing = @(x,t) 4*(3*x.^2 + 1) + 4*t*(3*r*x.^2 -(r-sigma^2/2)*6*x+(r-3*sigma^2));

bc_left = @(t) exact_sol(x_min,t);
bc_right = @(t) exact_sol(x_max,t);

initial_cond = @(x) 0;

% comparing solutions for Nt=512
Nt = 512;
[V,FD_grid,~] = bs_logform_forward_euler(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,Nx,T,Nt);

figure
set(gca,'FontSize',20)
plot(FD_grid,exact_sol(FD_grid,T),'ob','LineWidth',2,'DisplayName','Exact solution') 
hold on
plot(FD_grid,V(:,Nt+1),'r','LineWidth',2,'DisplayName','Finite diff. approx.')
 
grid on
legend show
set(legend,'Location','NorthWest')

% difference between exact solution and finite differences approximation
max(abs(V(:,Nt+1)'-exact_sol(FD_grid,T)))


% comparing solutions for Nt=256
Nt = 256;
[V,FD_grid,~] = bs_logform_forward_euler(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,Nx,T,Nt);

figure
set(gca,'FontSize',20)
plot(FD_grid,exact_sol(FD_grid,T),'ob','LineWidth',2,'DisplayName','Exact solution') 
hold on
plot(FD_grid,V(:,Nt+1),'r','LineWidth',2,'DisplayName','Finite diff. approx.')
 
grid on
legend show
set(legend,'Location','NorthWest')

% difference between exact solution and finite differences approximation
max(abs(V(:,Nt+1)'-exact_sol(FD_grid,T)))

%% Ex3

exact_sol = @(x,t) 4*t*(3*x.^2 + 1);

sigma = 1;
r = 1;
x_min = 0;
x_max = 1;
T = 1;
Nx = 21;
Nt = 100;
theta = 0.5;

forcing = @(x,t) 4*(3*x.^2 + 1) + 4*t*(3*r*x.^2 -(r-sigma^2/2)*6*x+(r-3*sigma^2));

bc_left = @(t) exact_sol(x_min,t);
bc_right = @(t) exact_sol(x_max,t);

initial_cond = @(x) 0;

[V,FD_grid,~] = bs_logform_timestepping(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,Nx,T,Nt,theta);

% comparing solutions
figure
set(gca,'FontSize',20)
plot(FD_grid,exact_sol(FD_grid,T),'ob','LineWidth',2,'DisplayName','Exact solution') 
hold on
plot(FD_grid,V(:,Nt+1),'r','LineWidth',2,'DisplayName','Finite diff. approx.')
 
grid on
legend show
set(legend,'Location','NorthWest')

% difference between exact solution and finite differences approximation
max(abs(V(:,Nt+1)'-exact_sol(FD_grid,T)))

%% Ex4

theta=0.5; % use crank-nicholson time-stepping scheme

T = 1;
sigma = 0.21;
r = 0.015;
K = 1;

Nx = 925;
Nt=500; 
deltat=T/Nt;

initial_cond = @(S) max(K-S,0);
V0 = @(x) initial_cond(exp(x));

forcing = @(x,t) 0*x;

tols = [1e-2 1e-3 1e-4 1e-5]; 
errors = [];

for tol=tols
    
    % compute the left and right boundary
    Smin = K*exp(-r*T -sigma^2*T/2+sigma*sqrt(T)*norminv(tol/K));
    Smax = K*exp(sigma^2*T/2 - sigma*sqrt(T)*norminv(tol/K));
    x_min = log(Smin);
    x_max = log(Smax);
    bc_left  = @(t) K*exp(-r*t)-Smin;
    bc_right = @(t) 0;
    
    % discretization parameters: observe that we need to keep h constant for every Smax
    h = (x_max-x_min)/Nx; 
    
    % finite differences solver
    [V,FD_grid,~] = bs_logform_timestepping(sigma,r,forcing,bc_left,bc_right,V0,x_min,x_max,Nx,T,Nt,theta);
    
     % exact solution with blsprice
    Vvect = exp(FD_grid); 
    Vexact = [];
    for v=Vvect
        [~,Put] = blsprice(v, K, r, T, sigma);
        Vexact = [Vexact Put];
    end
    
    % finite differences error
    errors = [errors max(abs(Vexact-V(:,end)'))];
    
    % evaluation of the finite differences solution for a given S value
    S_query = 0.9;
    V_query = interp1(FD_grid,V(:,end),log(S_query),'linear'); 
    [~,V_query_exact] = blsprice(S_query, K, r, T, sigma);
        
    % plot of the exact value and of the finite differences approximation
    % for tol = 1e-2
    if(tol==1e-2)
        figure
        set(gca,'FontSize',20) 
        plot(Vvect,Vexact,'ob','LineWidth',2,'DisplayName','Exact solution')
        hold on
        plot(Vvect,V(:,end),'r','LineWidth',2,'DisplayName','Finite Diff. approx.'); 
        grid on
        legend show
        set(legend,'Location','NorthEast')
        set(gca,'XTick',[Smin 1 1.5 Smax]) 
        set(gca,'XTickLabel',{'Smin','1','1.5','Smax'})
    end

end

figure
set(gca,'FontSize',20) 
loglog(tols,errors,'-xk','LineWidth',2,'MarkerSize',10,'DisplayName','truncation error') 
hold on
loglog(tols,tols,'--+','LineWidth',2,'MarkerSize',10,'DisplayName','tolerances')
grid on
legend show
set(legend,'Location','SouthEast')












