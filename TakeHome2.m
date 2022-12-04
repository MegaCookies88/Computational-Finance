% Data Loading :
load('Call_20050103.mat');
data = Call_20050103;

PCall    = data(:,1);       % call price
Strike   = data(:,2);       % strike price
Maturity = data(:,3)/252;   % maturity in years considering there are 252 trading days in average for the NYSE and NASDAQ
IV       = data(:,4);       % implied volatility


% Parameters Initialization :
r     = 0.015;
S     = 1202.10;
theta = 0.04;
kappa = 1.5;
sigma = 0.3;
rho   = -0.6;
V     = 0.0441;


% Optimization Step :

% optimization function 
d = @(x) VCall_Heston(Strike, Maturity, r, x(1), x(2), x(3), x(4), S, x(5)) - PCall; 
RMSE = @(x) sqrt(sum(d(x).^2)/length(PCall));

% initial values for optimization
x0 = [theta,kappa,sigma,rho,V]; 

% boundaries for variables
LB = [-inf,-inf,0,-1,0];
UB = [inf,inf,inf,0,inf]; % remark: we've seen in lecture 1 that there is empirical evidence that rho<0

% RMSE minimization
Calibrated_Parameters = fminsearchcon(RMSE,x0,LB,UB)
RMSE_min = RMSE(Calibrated_Parameters)


% Vectorized function of Call_Heston :
function [P] = VCall_Heston(K,T,r,nu,kappa,sigma,rho,S,V)
    P = ones(length(K),1);
    for i=1:length(K)
         P(i) = Call_Heston(K(i),T(i),r,nu,kappa,sigma,rho,S,V);
    end
end

