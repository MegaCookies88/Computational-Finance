function [V,FD_grid] = bs_logform_steady(sigma,r,rhs,bc_left,bc_right,x_min,x_max,Nx)

% [V,grid] = bs_logform_steady(sigma,r,rhs,bc_left,bc_right,x_min,x_max,Nx)
%
% solves the steady log-price version of the Black & Scholes model. 
%
% ----------------
% Inputs:
% ----------------
%
% sigma         : volatility (real number)
% r             : risk-free return (real number)
% rhs           : @-function describing the right-hand-side of the equation, rhs = @(x) ...
%                   x is a vector (spatial grid) and the output is a vector.
% bc_left       : value of the boundary condition on the left border (real value)
% bc_right      : value of the boundary condition on the right border (real value)
% x_min, x_max  : real values setting the extrema of the intervals in which the solution is computed
% Nx            : number of (space) intervals
%
% ----------------
% Outputs:
% ----------------
%
% V             : vector containing the solution of the PDE on each spatial node
% FD_grid       : spatial grid, contains the nodes on which the solution S is evaluated





% set grid and grid-size. 
h = (x_max-x_min)/Nx;
FD_grid = linspace(x_min,x_max,Nx+1); 
inner_grid = FD_grid(2:end-1);


% prepare rhs vector for the stationary solver, taking care of the contributions of the first and
% last node on the rhs vector coming from the discretization of transport (first derivative) and
% diffusion (second derivative)

ff = rhs(inner_grid);
ff(1) = ff(1) + 0.5*sigma^2*bc_left/h^2 - (r-sigma^2/2)*bc_left/(2*h); 
ff(end) = ff(end) + 0.5*sigma^2*bc_right/h^2 + (r-sigma^2/2)*bc_right/(2*h);


% prepare matrix for the stationary solver. We do not build the matrices but only their diagonals,
% since we will solve the system with thomas solver. Note the main diagonal
% has length Nx-1
% and the upper and lower have length Nx-2
% CHECK SLIDE 30 LECTURE 7

aa_up   = -0.5*sigma^2*ones(Nx-2,1)/h^2;
aa_main =  0.5*sigma^2*2*ones(Nx-1,1)/h^2;
aa_down = -0.5*sigma^2*ones(Nx-2,1)/h^2;

bb_up = -(r-sigma^2/2)*ones(Nx-2,1)/(2*h);
bb_down = (r-sigma^2/2)*ones(Nx-2,1)/(2*h);

cc_main = r*ones(Nx-1,1);

diag_up = aa_up+bb_up;
diag_main = aa_main+cc_main;
diag_down = aa_down+bb_down;


% solve and add boundary conditions. We use thomas solver
V = thomas(diag_down,diag_main,diag_up,ff);
V = [bc_left; V; bc_right];
