function [V,FD_grid,time_steps,iters_pSOR,errest_pSOR] = am_opt_pSOR(sigma,r,forcing,bc_left,bc_right,initial_cond,...
                                                                        S_min,S_max,Ns,T,deltat,theta,...
                                                                        tol,itermax,omega)


% [V,grid,time_steps] = am_opt_pSOR(sigma,r,forcing,bc_left,bc_right,initial_cond,S_min,S_max,Ns,T,deltat,theta)
%
% solves the american option pricing problem with a theta-method time-discretization scheme and a p-SOR
% algorithm for the solution of the linear complementarity problem
%
% ----------------
% Inputs:
% ----------------
%
% sigma         : @-function describing the volatility, sigma=@(S,t) ... It must accept vector values for S
% r             : risk-free return (real number)
% forcing       : @-function describing the right-hand-side of the equation, rhs = @(S,t) ...
%                   S is a vector (spatial grid), t is a real number and the output is a vector of the same dimension as S.
% bc_left       : @-function describing the boundary condition on the left border,
%                   bc_left = @(t) ...; t is a real number
% bc_right      : @-function describing the boundary condition on the right border,
%                   bc_right = @(t) ...; t is a real number
% initial_cond  : @-function describing the initial condition of the problem, initial_cond = @(S) ...
%                   S is a vector (spatial grid)    
% S_min, S_max  : real values setting the extrema of the intervals in which the solution is computed
% Ns            : the nodes of the spatial discretization are labeled 0, 1, ..., Ns, Ns+1. 
%                   Since Dirichlet boundary conditions on both sides are considered, we need to solve
%                   for Ns unknown (internal nodes)
% T             : final time of the equation
% deltat        : time-step
% theta         : parameter for the theta-method time-stepping scheme. 
%                   theta=0   --> Forward Euler
%                   theta=1   --> Backward Euler
%                   theta=0.5 --> Crank-Nicholson
%
%
% ----------------
% Outputs:
% ----------------
%
% V             : matrix containing the solution of the PDE. Each column represents the solution 
%                   on the whole spatial grid at a single time
% grid          : spatial grid, contains the nodes on which the solution S is evaluated
% time_steps    : time grid, contains the time steps at which S is computed



% set grid and grid-size.
h = (S_max-S_min)/(Ns+1);
FD_grid = linspace(S_min,S_max,Ns+2); 
inner_grid = FD_grid(2:end-1);

% set number of time steps
Nt=round(T/deltat);
time_steps = [0 deltat*(1:Nt)];
iters_pSOR=zeros(1,Nt);
errest_pSOR=zeros(1,Nt);

% init matrix containing the solution at each time step
V=zeros(length(FD_grid),Nt+1);
payoff=initial_cond(FD_grid)';
V(:,1)=payoff;


% the theta-method consists in the formula
%
% [ I + deltat*theta A ] V_new = V_old - deltat*(1-theta) A V_old  + deltat (theta *F_new +(1-theta) F_old ]
%
% we begin by building the matrix A, as in the steady case, and the matrix 
%
% B = [ I + deltat*theta A ] 
%
% while we need to build both A (as in the Forward Euler) and B (p-SOR)


% for reading convenience and for computational efficiency, we define some auxiliary vectors
S = inner_grid';
S_up = inner_grid(1:end-1)';
S_down = inner_grid(2:end)';

Ssq = inner_grid'.^2;
Ssq_up = inner_grid(1:end-1)'.^2;
Ssq_down = inner_grid(2:end)'.^2;



aa_up   = -0.5/h^2*sigma.^2.*ones(Ns-1,1).*Ssq_up; 
aa_main =  1/h^2*sigma.^2.*ones(Ns,1).*Ssq;
aa_down = -0.5/h^2*sigma.^2.*ones(Ns-1,1).*Ssq_down;

bb_up = -r*ones(Ns-1,1)/(2*h).*S_up;
bb_down = r*ones(Ns-1,1)/(2*h).*S_down;

cc_main = r*ones(Ns,1);

B_up = deltat*theta*(aa_up+bb_up);
B_main = deltat*theta*(aa_main+cc_main)+ones(Ns,1);
B_down = deltat*theta*(aa_down+bb_down);

A_up = aa_up+bb_up;
A_main = aa_main+cc_main;
A_down = aa_down+bb_down;


% A is zero outside the 3 central diagonals, therefore it's convenient to store it
% in sparse format. Note that we need to artificially extend the super-diagonal and sub-diagonal vectors
A = spdiags([NaN; A_up],1,Ns,Ns) + spdiags(A_main,0,Ns,Ns) + spdiags([A_down; NaN],-1,Ns,Ns);
B = spdiags([NaN; B_up],1,Ns,Ns) + spdiags(B_main,0,Ns,Ns) + spdiags([B_down; NaN],-1,Ns,Ns);


% time-stepping loop. We initialize the forcing with the corrections for the initial time-step
f_old = forcing(inner_grid,0)';
corr_left_old = 0.5*sigma^2*bc_left(0)/h^2*S(1).^2 - r*bc_left(0)/(2*h)*S(1);
corr_right_old = 0.5*sigma^2*bc_right(0)/h^2*S(end).^2 + r*bc_right(0)/(2*h)*S(end);
f_old(1) =  f_old(1) + corr_left_old;
f_old(end) =  f_old(end) + corr_right_old;

for tn = 1:Nt
    
    t_new=tn*deltat;

    % forcing at the new time-step
    f_new = forcing(inner_grid,t_new)';
    corr_left_new = 0.5*sigma^2*bc_left(t_new)/h^2*S(1).^2 - r*bc_left(t_new)/(2*h)*S(1);
    corr_right_new = 0.5*sigma^2*bc_right(t_new)/h^2*S(end).^2 + r*bc_right(t_new)/(2*h)*S(end);
    f_new(1) =  f_new(1) + corr_left_new;
    f_new(end) =  f_new(end) + corr_right_new;

    % we now need to solve the LCP
    %
    % BU-f_tot>=0
    % U-payoff>=0
    % (BU-f_tot)'(U-payoff)=0
    %
    % assemble f_tot and solve the LCP using pSOR method
    %
    % f_tot = V_old - deltat*(1-theta) A V_old  + deltat (theta *F_new +(1-theta) F_old ]
   
    f_tot = V(2:end-1,tn) - deltat*(1-theta)*A*V(2:end-1,tn) + deltat*theta*f_new + deltat*(1-theta)*f_old;
        
    % solve with pSOR, using the solution at the previous time-step as initial guess
    [sol,k,errest]=pSOR(B,f_tot,omega,payoff(2:end-1),tol,itermax,V(2:end-1,tn));
    iters_pSOR(tn)=k;    
    errest_pSOR(tn)=errest;
    
    V(:,tn+1)=[bc_left(t_new); sol; bc_right(t_new)]; 
    % update the forcing at the previous time step
    f_old = f_new;
end
