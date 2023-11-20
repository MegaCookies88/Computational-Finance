function [I_N] = ChebInterpol(f,x,n,a,b)
% Computes the Chebyshev interpolation of order n of the function f

% Input:
% f: function on vector
% x: points 
% n: order
% [a,b]: interval (optional)

% Output:
% I: interpolated prices

if ~exist('a','var')
    a = min(x);
end

if ~exist('b','var')
    b = max(x);
end

% Initializing
I_N = zeros(1,length(x));

% Scaling
x = (x-a)./(b-a)*2-1;

% Nodes for [-1,1] then scaling for [a,b]
p = cos(pi*(0:n)/n);
p = 0.5*(a+b)+0.5*(b-a)*p;

% Weights
c = zeros(1,n+1);

try 
    price = f(p);
catch
    price = zeros(size(p));
    for i=1:size(p,2)
        price(i) = f(p(i));
    end
end

for i=0:n
    c(i+1) = 2.^(i>0)/n * dot([0.5 ones(1,n-1) 0.5].*price, cos(i*pi*(0:n)/n));
end

% Interpolated price
for i=1:length(x)
    T = cos((0:n)*acos(x(i)));
    I_N(i) = dot(c,T);
end

end

