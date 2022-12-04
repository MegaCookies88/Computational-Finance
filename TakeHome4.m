%%% Part b %%%

% PDE parameters
sigma = @(S,t) 0.2;
r = 0.015;
K = 1;
T = 1;
S_max = 2.5; 

% theta-method parameters
Nt = 100;
Nx = 60;
theta = 0.5;
forcing = @(S,t) 0;
initial_cond = @(S) max(0,K-S);
bc_right = @(t) 0;


%%% Question i  %%%

% Numerical solution for PDE (3)
[V,FD_grid,time_steps] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);

% -------------------------------------------------------------------------

%%% Question ii  %%%

% Numerial solution for PDE of Vega
h = S_max/Nx; % step size for S
dt = T/Nt; % step size for T

forcing = @(S,t) sigma(S,t) * S.^2 .* V_SS(S,t,h,V,S_max,dt);
initial_cond = @(S) 0;

[Vega,FD_grid_vega,time_step_vega]= bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);

% -------------------------------------------------------------------------

%%% Question iii  %%%

% Comparing with exact solution in a plot
exact_vega = zeros(1,Nx);
for i=2:Nx
       exact_vega(i)= blsvega(FD_grid_vega(i),K,r,1,sigma(1,1));
end

figure(1)
plot(FD_grid_vega(2:end), Vega(1:end-1,end))
hold on
plot(FD_grid_vega(2:end), exact_vega)
xlabel('S')
ylabel('\nu')
title('Comparison of the exact and the approximated \nu')
legend("approximate \nu","exact \nu")
hold off
 
% -------------------------------------------------------------------------

%%% Part c %%%

h_tab = zeros(1,6);
error = zeros(1,6);

for i=0:5
    
    Nx = 10*2^i;
    Nt = 0.6*Nx;
    
    h = S_max/Nx;
    dt = T/Nt;
    
    h_tab(i+1) = h;
    
    % Numerical solution for PDE (3)
    forcing = @(S,t) 0;
    initial_cond = @(S) max(0,K-S);
    bc_right = @(t) 0;
    
    [V,FD_grid,time_step] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);
    
    % Numerial solution for PDE of Vega
    forcing = @(S,t) sigma(S,t) * S.^2 .* V_SS(S,t,h,V,S_max,dt);
    initial_cond = @(S) 0;

    [Vega,FD_grid_vega,time_step_vega]= bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);
    
    % Computing exact solution
    exact_vega = zeros(1,Nx);
    for j=2:Nx
           exact_vega(j)= blsvega(FD_grid_vega(j),K,r,1,sigma(1,1));
    end
    
    error(i+1)= max(abs(Vega(1:end-1,end)'-exact_vega));

end


% Plotting error 
figure(2)
plot(h_tab,error,'r--')
hold on
xlabel('h')
ylabel('error')
title('Comparison of errors of approximation')

% -------------------------------------------------------------------------

%%% Part d %%%

% PDE parameters
sigma = @(S,t) 0.2;
r = 0.015;
K = 1;
T = 1;
S_max = 2.5; 

% theta-method parameters
Nt = 100;
Nx = 60;
theta = 0.5;
forcing = @(S,t) 0;
initial_cond = @(S) max(0,K-S);
bc_right = @(t) 0;


%%% Question i  %%%

d_s_tab = [0.0005 0.001 0.01 0.02 0.03 0.05];
vega_tab = {}; % cell that saves Vega for each delta sigma in the tab

% Computing new approximations of Vega
for i=1:6

    d = d_s_tab(i);
    
    sigma = @(S,t) 0.2-d;
    [Vl] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);
    
    sigma = @(S,t) 0.2+d;
    [Vr] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);

    vega_tab(i,:) = {d,(Vr-Vl)/(2*d)};

end

% -------------------------------------------------------------------------

%%% Question ii  %%%

% Computing the errors
h_tab = zeros(1,6);
error = zeros(6,6);

for i=0:5
    
    Nx = 10*2^i;
    Nt = 0.6*Nx;
    
    h = S_max/Nx;
    dt = T/Nt;
    
    h_tab(i+1) = h;

    % Computing error for each delta-sigma
    for j=1:6
        
        % Computing approximation
        d = d_s_tab(j);
        sigma = @(S,t) 0.2-d;
        [Vl] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);
        sigma = @(S,t) 0.2+d;
        [Vr,FD_grid_vega] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta);
        
        Vega = (Vr-Vl)/(2*d);
        
        % Computing exact solution
        sigma = @(S,t) 0.2;
        exact_vega = zeros(1,Nx);
        for l=2:Nx
           exact_vega(l)= blsvega(FD_grid_vega(l),K,r,1,sigma(1,1));
        end
        
        % Error
        error(j,i+1)= max(abs(Vega(1:end-1,end)'-exact_vega));
        
    end

end


% Plotting errors on the previous graph
for j=1:6
    plot(h_tab,error(j,:))
end
fplot(@(x) (0.3*x),[0,0.25],'b--')
legend("V_{SS} approximation (b)",...
       "\delta\sigma=0.0005","\delta\sigma=0.001",...
       "\delta\sigma=0.01","\delta\sigma=0.02",...
       "\delta\sigma=0.03","\delta\sigma=0.05",...
       "0.3h",...
       'Location','northwest')
xlabel('h')
ylabel('error')
title('Comparison of errors of approximation')
hold off


function vss = V_SS(S,t,h,V,S_max,dt)   
% this function computes the second derivative of V at (S,t)

vss = zeros(1,length(S));

for i=1:length(S)
    
if(S(i)<S_max && S(i)>0)
          vss(i) = (V(floor(S(i)./h),floor(t/dt)+1)...
                        -2*V(floor(S(i)./h)+1,floor(t/dt)+1)...
                        +V(floor(S(i)./h)+2,floor(t/dt)+1))...
                        /h^2;
else
    vss(i) = 0;
end

end

end




