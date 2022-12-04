function [V,FD_grid,time_steps] = bs_logform_timestepping(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,Nx,T,Nt,theta)


% [V,grid,time_steps] = bs_logform_timestepping(sigma,r,forcing,bc_left,bc_right,initial_cond,x_min,x_max,N_tot,T,deltat,theta)
%
% solves the time-dependent log-price version of the Black & Scholes model. 
%
% ----------------
% Inputs:
% ----------------
%
% sigma         : volatility (real number)
% r             : risk-free return (real number)
% forcing       : @-function describing the right-hand-side of the equation, rhs = @(x,t) ...
%                   x is a vector (spatial grid), t is a real number and the output is a vector of the same dimension as x.
% bc_left       : @-function describing the boundary condition on the left border,
%                   bc_left = @(t) ...; t is a real number
% bc_right      : @-function describing the boundary condition on the right border,
%                   bc_right = @(t) ...; t is a real number
% initial_cond  : @-function describing the initial condition of the problem, initial_cond = @(x) ...
%                   x is a vector (spatial grid)    
% x_min, x_max  : real values setting the extrema of the intervals in which the solution is computed
% Nx            : nuber of space intervals
% T             : final time of the equation
% Nt            : number of time intervals
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
h = (x_max-x_min)/Nx;
FD_grid = linspace(x_min,x_max,Nx+1); 
inner_grid = FD_grid(2:end-1);


% set time step
deltat=T/Nt;
time_steps = [0 deltat*(1:Nt)];

% init matrix containing the solution at each time step
V=zeros(length(FD_grid),Nt+1);
V(:,1)=initial_cond(FD_grid)';


% the theta-method consists in the formula
%
% [ I + deltat*theta A ] V_new = V_old - deltat*(1-theta) A V_old  + deltat (theta *F_new +(1-theta) F_old ]
%
% we begin by building the matrix A, as in the steady case, and the matrix 
%
% B = [ I + deltat*theta A ] 
%
% We will solve the linear system in B, hence we only need its diagonals, while we need
% to build A (as in the Forward Euler)


aa_up   = -0.5*sigma^2*ones(Nx-2,1)/h^2;
aa_main =  0.5*sigma^2*2*ones(Nx-1,1)/h^2;
aa_down = -0.5*sigma^2*ones(Nx-2,1)/h^2;

bb_up = -(r-sigma^2/2)*ones(Nx-2,1)/(2*h);
bb_down = (r-sigma^2/2)*ones(Nx-2,1)/(2*h);

cc_main = r*ones(Nx-1,1);

B_up = deltat*theta*(aa_up+bb_up);
B_main = deltat*theta*(aa_main+cc_main)+ones(Nx-1,1);
B_down = deltat*theta*(aa_down+bb_down);

A_up = aa_up+bb_up;
A_main = aa_main+cc_main;
A_down = aa_down+bb_down;



% A is zero outside the 3 central diagonals, therefore it's convenient to store it
% in sparse format. Note that we need to artificially extend the super-diagonal and sub-diagonal vectors
A = spdiags([NaN; A_up],1,Nx-1,Nx-1) + spdiags(A_main,0,Nx-1,Nx-1) + spdiags([A_down; NaN],-1,Nx-1,Nx-1);


% time-stepping loop. We initialize the forcing with the corrections for the initial time-step
f_old = forcing(inner_grid,0)';
corr_left_old = 0.5*sigma^2*bc_left(0)/h^2 - (r-sigma^2/2)*bc_left(0)/(2*h);
corr_right_old = 0.5*sigma^2*bc_right(0)/h^2 + (r-sigma^2/2)*bc_right(0)/(2*h);
f_old(1) =  f_old(1) + corr_left_old;
f_old(end) =  f_old(end) + corr_right_old;

for tn = 1:Nt
    
    t_new=tn*deltat;

    % forcing at the new time-step
    f_new = forcing(inner_grid,t_new)';
    corr_left_new = 0.5*sigma^2*bc_left(t_new)/h^2 - (r-sigma^2/2)*bc_left(t_new)/(2*h);
    corr_right_new = 0.5*sigma^2*bc_right(t_new)/h^2 + (r-sigma^2/2)*bc_right(t_new)/(2*h);
    f_new(1) =  f_new(1) + corr_left_new;
    f_new(end) =  f_new(end) + corr_right_new;


    % we now need to solve the linear system
    %
    % (I + deltat*theta A) V_new = f_tot
    %
    % assemble f_tot and solver the linear system using Thomas
    %
    % f_tot = V_old - deltat*(1-theta) A V_old  + deltat (theta *F_new +(1-theta) F_old ]

    f_tot = V(2:end-1,tn) - deltat*(1-theta)*A*V(2:end-1,tn) + deltat*theta*f_new + deltat*(1-theta)*f_old;
        
    % solve the system and store the solution
    sol = thomas(B_down,B_main,B_up,f_tot);
    V(:,tn+1)=[bc_left(t_new); sol; bc_right(t_new)];

    % update the forcing at the previous time step
    f_old = f_new;
end
