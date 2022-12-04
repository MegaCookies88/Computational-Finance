%% EX1

D = 2 * eye(5) - diag(ones(4,1),1) - diag(ones(4,1),-1);

x = randn(10,1);
y = randn(10,1);
z = x.*y;
x'*y - sum(z);

A = eye(4);
A = A(:,[4,3,2,1]);

func(x,y);

A = rand(1024); B = rand(1024); C = zeros(1024); D = zeros(1024);
disp('Time needed for the built−in matrix multiplication'); tic; C = A*B; toc; 
%disp('Time needed for the function matmult'); tic; D = matmult(A,B); toc;

f = @(x) exp(-x.*x);
quad(f,-5,5)^2 

p = @(x) x.^4+10*x.^3-372*x.^2-2714*x+20995
x = -20:1/100:20;
r = roots([1 10 -372 -2714 20995]);
plot(x,p(x))
hold on
plot(r,zeros(1,size(r,1)),'o')
hold off


%% EX2

f = @(x) 3./(x+1);
g = @(x) x+3;

x = linspace(-10,10,100);
hold on
plot(x,f(x),'Red')
plot(x,g(x),'Blue')
legend('f(x)', 'g(x)')
hold off

h = @(x) f(x)-g(x);
bisection(h,-5,2,1e-3)


%% EX3

figure(1);
x = linspace(0.1, 100, 1000);
plot(x, f(x), 'b', [0 100], [0 0], 'r');
figure(2);
x = linspace(0.5, 10, 1000);
plot(x, f(x), 'b', [0 10], [0 0], 'r');

for i=1:10
    fprintf('start from %d\n', i);
    newton(@f, @df, i, 1e-6); 
end


%% EX4
  
SN = @(N,a,b,x) 0.5*a(0) + a(2:N+1)'*cos((1:1:N)*x) + b(2:N+1)'*sin((1:1:N)*x); 
x = linspace(-pi,pi,50);
%hold on
plot(x,ff(x))


%% Functions always declared last
function C = matmult(A,B)
%multiply A and B
C = zeros(size(A,1),size(B,2));
for i = 1:size(A,1)
    for j = 1:size(B,2)
        C(i,j) = 0;
        for l = 1:size(A,2)
            C(i,j) = C(i,j) + A(i,l) * B(l,j);
        end
    end
end
end

function x = newton(f,df,x0,tol)
% Newton method for computing a zero of a function f. % Input:
% f − function handle to f
% df − function handle to f'
% x − starting value
% tol − User tolerance
% Output:
% x0 − approximate zero of f
it = 0; maxit = 100;
x = x0;
while (it < maxit & abs(f(x)) >= tol)
    it = it + 1;
    x = x - f(x) / df(x);
    fprintf(' step %d: x=%f, |f(x)|=%e\n', it-1, x, abs(f(x)));
end
if (abs(f(x)) >= tol)
    error('Convergence failure.');
end
end

function y = f(x) 
if x>0
    y = 2./(exp(x)-1) -2./sqrt(x)+1;
else
    y = 0;
end
end

function y = df(x)
if x > 0
    y = x.^(-1.5) - 2*exp(x)./(exp(x)-1).^2; 
else
    y = 0;
end
end

function y = ff(x)
if x >= -pi & x <= 0    
    y = 1;
end
if x <= pi & x > 0
    y = -1;
end 
y = 0;
end

function a = Coeffa(f,N)
a = zeros(1,N+1);
a(1) = 1/pi * quad(f,-pi,pi);
for n = 1:1:N
    fn = @(x) cos(n * x) .* f(x);
    a(n+1) = 1/pi * quad(fn,-pi,pi);
end
end

function b = Coeffb(f,N)
b = zeros(1,N+1);
for n = 1:1:N
    fn = @(x) sin(n * x) .* f(x);
    b(n+1) = 1/pi * quad(fn,-pi,pi);
end
end


