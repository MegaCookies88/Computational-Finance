function zero = bisection(f,a,b,tol)

% Input: function handle f, interval limits a and b, tolerance
% Output: the algorithm gives a solution to the equation f(x)=0
%         on the interval [a,b] within the given tolerance by 
%         using the bisection method


if f(a)*f(b)>0 
    disp('Initial condition not satisfied')
else
    x = (a + b)/2;
    err = abs(f(x));
    while err > tol
        if f(a)*f(x)<0 
             b = x;
        else
             a = x;          
        end
    x = (a + b)/2;
    err = abs(f(x));
    fprintf('x= %f, |f(x)|=%e\n',x, err)
    end
    zero = x;
end