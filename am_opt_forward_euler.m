function [V,FD_grid,time_steps] = am_opt_forward_euler(sigma,r,forcing,bc_left,bc_right,initial_cond,S_min,S_max,Ns,T,deltat)


% [V,grid,time_steps] = am_opt_forward_euler(sigma,r,forcing,bc_left,bc_right,initial_cond,S_min,S_max,Ns,T,deltat)
%
% solves the american option pricing problem with a Forward Euler time-discretization scheme
%
% ----------------
% Inputs:
% ----------------
%
% sigma         : volatility (real number)
% r             : risk-free return (real number)
% forcing       : @-function describing the right-hand-side of the equation, forcing = @(S,t) ...
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
%
%
% ----------------
% Outputs:
% ----------------
%
% V             : matrix containing the solution of the PDE. Each column represents the solution 
%                   on the whole spatial grid at a single time
% FD_grid       : spatial grid, contains the nodes on which the solution S is evaluated
% time_steps    : time grid, contains the time steps at which S is computed



% set grid and grid-size. 
h = (S_max-S_min)/(Ns+1);
FD_grid = linspace(S_min,S_max,Ns+2); 
inner_grid = FD_grid(2:end-1);


% set number of time steps
Nt=round(T/deltat);
time_steps = [0 deltat*(1:Nt)];

% init matrix containing the solution at each time step
V=zeros(length(FD_grid),Nt+1);
payoff=initial_cond(FD_grid)';
V(:,1)=payoff;


% the american option Forward Euler method consists in the formula
%
% V_t = V_old - deltat A V_old  + deltat F_old
% V_new = max(V_t, init_cond(S) )
%
% we begin by building the matrix A, as in the steady case.

% for reading convenience and for computational efficiency, we define some auxiliary vectors
S = inner_grid';
S_up = inner_grid(1:end-1)';
S_down = inner_grid(2:end)';

Ssq = inner_grid'.^2;
Ssq_up = inner_grid(1:end-1)'.^2;
Ssq_down = inner_grid(2:end)'.^2;


aa_up   = -0.5*sigma^2*ones(Ns-1,1)/h^2.*Ssq_up;
aa_main =  0.5*sigma^2*2*ones(Ns,1)/h^2.*Ssq;
aa_down = -0.5*sigma^2*ones(Ns-1,1)/h^2.*Ssq_down;

bb_up = -r*ones(Ns-1,1)/(2*h).*S_up;
bb_down = r*ones(Ns-1,1)/(2*h).*S_down;

cc_main = r*ones(Ns,1);

diag_up = aa_up+bb_up;
diag_main = aa_main+cc_main;
diag_down = aa_down+bb_down;

% A consists mainly of zeros, so we store it in sparse format. Note that
% we need to artificially extend the subdiagonal and superdiagonal vectors:
% add 1 value at te beginning of the super-diagonal and 1 value at the end of the sub-diagonal
A = spdiags([NaN; diag_up],1,Ns,Ns) + spdiags(diag_main,0,Ns,Ns) + spdiags([diag_down; NaN],-1,Ns,Ns);


% time-stepping loop
t_old=0;

for tn = 1:Nt
    
    t_new=t_old+deltat;
    
    % build update vector according to the formula
    %
    % V_new = V_old - deltat A V_old  + deltat F_old
    %
    % note that V_old is already stored in a suitable column of the matrix V

    % we first build F_old with boundary corrections (evaluated ad t_old, according to Forward Euler formula)
    F_old = forcing(inner_grid,t_old)';
    corr_left = 0.5*sigma^2*bc_left(t_old)/h^2*S(1).^2 - r*bc_left(t_old)/(2*h)*S(1);
    F_old(1) = F_old(1)+corr_left;
    corr_right = 0.5*sigma^2*bc_right(t_old)/h^2*S(end).^2 + r*bc_right(t_old)/(2*h)*S(end);
    F_old(end) = F_old(end)+corr_right;
    
    % then update solution according to the formula above 
    Vnew = V(2:end-1,tn) - deltat*A*V(2:end-1,tn) + deltat*F_old;
    
    % apply the constraint and store the solution
    V(:,tn+1)=max([bc_left(t_new); Vnew; bc_right(t_new)],payoff);
    
    t_old=t_new;
end
