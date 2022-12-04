Strikes  = [50;80;90;95;100;105;110;120;150];
Maturity = [1/12;1/6;1/4;1/2;1;2;3;5];
[T,K] = meshgrid(Maturity, Strikes);

S       = 100;
r       = 0.1;
theta   = -0.14;
nu      = 0.2;
sigma   = 0.12;

alpha   = -2; % if alpha < -1 gives put prices using the call formula

nMat = length(Maturity); 
nStrikes = length(Strikes);
d = length(K(:));

P  = zeros(d,1); % matrix of prices using FFT scheme
IV = zeros(d,1); % matrix of implied volatilites

% Using trapezoid rule
for i=1:d
    P(i) = Call_FFT_VG(K(i), T(i), r, theta, nu, sigma, S, alpha, 0);
    
    if P(i)<=0
        P(i) = eps;
    end
    
    if alpha>0 % call options
        IV(i) = blsimpv(S, K(i), r, T(i), P(i),[], 0, [], true);
    end
    
    if alpha<0 % put options
        IV(i) = blsimpv(S, K(i), r, T(i), P(i),[], 0, [], false);
    end
    
end

IV = reshape(IV, nStrikes, nMat);
P = reshape(P, nStrikes, nMat);
mesh(IV);
axis tight;
set(gca,'XTickLabel',num2str(Maturity));
set(gca,'YTickLabel',num2str(Strikes));
xlabel('Maturities');
ylabel('Strikes');
title('Implied Vol surface');

disp(P)

% Using simpson rule
for i=1:d
    P(i) = Call_FFT_VG(K(i), T(i), r, theta, nu, sigma, S, alpha, 1);
    
    if P(i)<=0
        P(i) = eps;
    end
    
    if alpha>0 % call options
        IV(i) = blsimpv(S, K(i), r, T(i), P(i),[], 0, [], true);
    end
    
    if alpha<0 % put options
        IV(i) = blsimpv(S, K(i), r, T(i), P(i),[], 0, [], false);
    end
    
end

IV = reshape(IV, nStrikes, nMat);
P = reshape(P, nStrikes, nMat);
mesh(IV);
axis tight;
set(gca,'XTickLabel',num2str(Maturity));
set(gca,'YTickLabel',num2str(Strikes));
xlabel('Maturities');
ylabel('Strikes');
title('Implied Vol surface');

disp(P)
