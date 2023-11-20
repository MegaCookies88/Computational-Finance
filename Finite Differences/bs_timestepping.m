function [V,FD_grid,time_steps] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond,S_max,Nx,T,Nt,theta)

% [V,grid,time_steps] = bs_timestepping(sigma,r,forcing,bc_right,initial_cond S_max,Nx,T,Nt,theta)
%
% solves the time-dependent version of the Black & Scholes model in the interval [0,Smax]. This meaNx that no boundary conditioNx
% on 0 need to be prescribed since in 0 the equation is  \partial_t u + rU=0, U(0)=initial_cond(0) i.e. we are integrating an
% ODE. We modify the meaning of the input Nx, which is now intended as the number of intervals for the
% discretization of [0,Smax]. This meaNx that we have Nx+1 nodes and we solve for Nx unknowNx since the node
% corresponding to Smax has a Dirichlet Boundary Condition.
%
% ----------------
% Inputs:
% ----------------
%
% sigma         : @-function describing the volatility, sigma=@(S,t) ... It must accept vector values for S
% r             : risk-free return (real number)
% forcing       : @-function describing the right-hand-side of the equation, rhs = @(S,t) ...
%                   S is a vector (spatial grid), t is a real number and the output is a vector of the same dimeNxion as S.
% bc_left       : @-function describing the boundary condition on the left border,
%                   bc_left = @(t) ...; t is a real number
% bc_right      : @-function describing the boundary condition on the right border,
%                   bc_right = @(t) ...; t is a real number
% initial_cond  : @-function describing the initial condition of the problem, initial_cond = @(S) ...
%                   S is a vector (spatial grid)    
% S_min, S_max  : real values setting the extrema of the intervals in which the solution is computed
% Nx            : number of intervals in the discretization of [0,Smax]
% T             : final time of the equation
% Nt            : number of time-steps
% theta         : parameter for the theta-method time-stepping scheme. 
%                   theta=0   --> Forward Euler
%                   theta=1   --> Backward Euler
%                   theta=0.5 --> Crank-Nicholson
%
% ----------------
% Outputs:
% ----------------
%
% V             : matrix containing the solution of the PDE. Each column represents the solution 
%                   on the whole spatial grid at a single time
% grid          : spatial grid, contaiNx the nodes on which the solution S is evaluated
% time_steps    : time grid, containing the time steps at which S is computed


S_min=0;

% set grid and grid-size.
h = (S_max-S_min)/Nx;
FD_grid = linspace(S_min,S_max,Nx+1); 
inner_grid = FD_grid(1:end-1);

% set number of time steps
deltat=T/Nt;
time_steps = [0 deltat*(1:Nt)];

% init matrix containing the solution at each time step
V=zeros(length(FD_grid),Nt+1);
V(:,1)=initial_cond(FD_grid)';

% for reading convenience and for computational efficiency, we define some auxiliary vectors
S = inner_grid';
S_up = inner_grid(1:end-1)'; 
S_down = inner_grid(2:end)'; 

Ssq = inner_grid'.^2;
Ssq_up = inner_grid(1:end-1)'.^2;
Ssq_down = inner_grid(2:end)'.^2;


% since sigma depends now on t, we need to concider a theta-method with time-dependent A, 
% i.e. the following formula
%
% B V_new = V_old - deltat*(1-theta) A_old V_old  + deltat (theta*F_new +(1-theta)*F_old),
% B = [ I + deltat*theta A_new ] 
%
% note also that now B, A_old contain a term depending on S. However, we can build some
% parts of A and B before the temporal loop, for computational efficiency.


aa_old_up   = -0.5/h^2*sigma(S_up,0).^2.*Ssq_up; 
aa_old_main =  1/h^2*sigma(S,0).^2.*Ssq;
aa_old_down = -0.5/h^2*sigma(S_down,0).^2.*Ssq_down;

bb_up = -r*ones(Nx-1,1)/(2*h).*S_up;
bb_down = r*ones(Nx-1,1)/(2*h).*S_down;

cc_main = r*ones(Nx,1);


% time-stepping loop. We initialize the forcing with the correctioNx for the initial time-step
f_old = forcing(inner_grid,0)';
corr_right_old = 0.5*sigma(S(end),0)^2*bc_right(0)/h^2*S(end).^2 + r*bc_right(0)/(2*h)*S(end);
f_old(end) =  f_old(end) + corr_right_old;

for tn = 1:Nt
    
    t_new=tn*deltat;

    % forcing at the new time-step
    f_new = forcing(inner_grid,t_new)';
    corr_right_new = 0.5*sigma(S(end),t_new)^2*bc_right(t_new)/h^2*S(end).^2 + r*bc_right(t_new)/(2*h)*S(end);
    f_new(end) =  f_new(end) + corr_right_new;

    
    % we now need to solve the linear system
    %
    % (I + deltat*theta A) V_new = f_tot
    %
    % with
    %
    % f_tot = V_old - deltat*(1-theta) A V_old  + deltat (theta*f_new +(1-theta)*f_old)
   

    % A_old is the matrix multiplying V_old on the rhs of the linear system. It depends on t_old, hence it changes
    % at every time-step, therefore it is not possible to build it before the time loop
    A_old_up   = aa_old_up+bb_up;
    A_old_main = aa_old_main+cc_main;
    A_old_down = aa_old_down+bb_down;
    A_old = spdiags([NaN; A_old_up],1,Nx,Nx) + spdiags(A_old_main,0,Nx,Nx) + spdiags([A_old_down; NaN],-1,Nx,Nx);

    f_tot = V(1:end-1,tn) - deltat*(1-theta)*A_old*V(1:end-1,tn) + deltat*theta*f_new + deltat*(1-theta)*f_old;

    
    % now build the diagonals for A new and B
    aa_new_up   = -0.5/h^2*sigma(S_up,t_new).^2.*Ssq_up;
    aa_new_main =  1/h^2*sigma(S,t_new).^2.*Ssq;
    aa_new_down = -0.5/h^2*sigma(S_down,t_new).^2.*Ssq_down;
        
    B_up = deltat*theta*(aa_new_up+bb_up);
    B_main = deltat*theta*(aa_new_main+cc_main)+ones(Nx,1);
    B_down = deltat*theta*(aa_new_down+bb_down);
        
         
    % solve the system and store the solution
    sol = thomas(B_down,B_main,B_up,f_tot);
    V(:,tn+1)=[sol; bc_right(t_new)];

    % update the forcing at the previous time step
    f_old = f_new;
    
    % update aa_* vectors
    aa_old_up = aa_new_up;
    aa_old_main = aa_new_main;
    aa_old_down = aa_new_down;

end
