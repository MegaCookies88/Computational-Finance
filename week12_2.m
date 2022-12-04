X = linspace(-5,5,100);
mu = zeros(size(X));
num_samp = 5;

% squared exponential kernel
theta=[1,0.3];
sigma=sqrdexp(X',X',theta(1),theta(2));
plot_GP(X,mu,sigma,num_samp)

% linear kernel
theta=[1,0.3,0.5];
sigma=linear_kernel(X',X',theta(1),theta(2),theta(3)); 
plot_GP(X,mu,sigma,num_samp)

% periodic kernel
theta=[1,3,pi];
sigma=periodic_kernel(X',X',theta(1),theta(2),theta(3));
plot_GP(X,mu,sigma,num_samp)

% parameter effect
for i=linspace(0.1,1,3)
    theta=[1,i];
    sigma=sqrdexp(X',X',theta(1),theta(2));
    %plot_GP(X,mu,sigma,num_samp)
end

%--------------------------------------------------------------------------

% Kernels functions

% Squared Exponential
function [cov] = sqrdexp(Xn,Xm,sigma,length_scale)
    l = length_scale;
    cov = (sigma^2)*exp(-(pdist2(Xn,Xm).^2)/(2*l^2));
end

% Linear
function [cov] = linear_kernel(Xn,Xm,sigma0,sigma1,c)
    cov = (sigma0^2)+(sigma1^2)*(Xn-c)*(Xm-c)';
end

% Periodic
function [cov] = periodic_kernel(Xn,Xm,sigma,length_scale,p)
    l = length_scale;
    cov = (sigma^2)*exp(-(2/(l^2))*(sin(pi*pdist2(Xn,Xm)/p).^2));
end

% Plots the samples from the Gaussian Prior
function [] = plot_GP(X,mu,sigma,num_samp)
figure()
    Y1 = mvnrnd(mu,sigma,num_samp);
    for i =1:num_samp 
        plot(X,Y1(i,:))
        hold on
    end
    plot(X,mu')
    title('samples functions')
    hold off
end