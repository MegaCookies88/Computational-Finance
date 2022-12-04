function [V,FD_grid,time_steps] = bs_logform_forward_euler(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,Nx,T,Nt)


% [V,grid,time_steps] = bs_logform_forward_euler(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,Nx,T,deltat)
%
% solves the time-dependent log-price version of the Black & Scholes model with Forward Euler
%
% ----------------
% Inputs:
% ----------------
%
% sigma         : volatility (real number)
% r             : risk-free return (real number)
% forcing       : @-function describing the right-hand-side of the equation, forcing = @(x,t) ...
%                   x is a vector (spatial grid), t is a real number and the output is a vector of the same dimension as x.
% bc_left       : @-function describing the boundary condition on the left border,
%                   bc_left = @(t) ...; t is a real number
% bc_right      : @-function describing the boundary condition on the right border,
%                   bc_right = @(t) ...; t is a real number
% initial_cond  : @-function describing the initial condition of the problem, initial_cond = @(x) ...
%                   x is a vector (spatial grid)    
% x_min, x_max  : real values setting the extrema of the intervals in which the solution is computed
% Nx            : number of space intervals
% Nt            : number of time intervals
% T             : final time
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
h = (x_max-x_min)/(Nx);
FD_grid = linspace(x_min,x_max,Nx+1); 
inner_grid = FD_grid(2:end-1);


% set number of time steps
deltat=T./Nt;
time_steps = [0 deltat*(1:Nt)];

% init matrix containing the solution at each time step
V=zeros(length(FD_grid),Nt+1);
V(:,1)=initial_cond(FD_grid)';


% the Forward Euler method consists in the formula
%
% V_new = V_old - deltat A V_old  + deltat F_old
%
% we begin by building the matrix A, as in the steady case.


aa_up   = -0.5*sigma^2*ones(Nx-2,1)/h^2;
aa_main =  0.5*sigma^2*2*ones(Nx-1,1)/h^2;
aa_down = -0.5*sigma^2*ones(Nx-2,1)/h^2;

bb_up = -(r-sigma^2/2)*ones(Nx-2,1)/(2*h);
bb_down = (r-sigma^2/2)*ones(Nx-2,1)/(2*h);

cc_main = r*ones(Nx-1,1);

diag_up = aa_up+bb_up;
diag_main = aa_main+cc_main;
diag_down = aa_down+bb_down;

% A consists mainly of zeros, so we store it in sparse format. Note that
% we need to artificially extend the subdiagonal and superdiagonal vectors:
% add 1 value at te beginning of the super-diagonal and 1 value at the end of the sub-diagonal
A = spdiags([NaN; diag_up],1,Nx-1,Nx-1) + spdiags(diag_main,0,Nx-1,Nx-1) + spdiags([diag_down; NaN],-1,Nx-1,Nx-1);


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
    corr_left = 0.5*sigma^2*bc_left(t_old)/h^2 - (r-sigma^2/2)*bc_left(t_old)/(2*h);
    F_old(1) = F_old(1)+corr_left;
    corr_right = 0.5*sigma^2*bc_right(t_old)/h^2 + (r-sigma^2/2)*bc_right(t_old)/(2*h);
    F_old(end) = F_old(end)+corr_right;
    
    % then update solution according to the formula above 
    Vnew = V(2:end-1,tn) - deltat*A*V(2:end-1,tn) + deltat*F_old;
    
    % store the solution
    V(:,tn+1)=[bc_left(t_new); Vnew; bc_right(t_new)];
    
    t_old=t_new;
end
