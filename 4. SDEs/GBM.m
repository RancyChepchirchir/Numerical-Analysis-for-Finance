% Simulate a Geometric Brownian motion
% dX = mu*X*dt + sigma*X*dW

% Define parameters and time grid
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
mu = 0.2; sigma = 0.3; % model parameters

% Compute the increments of an arithmetic Brownian motion
dlogX = (mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*randn(nsteps,npaths);

% Accumulate the increments
logX = [zeros(1,npaths); cumsum(dlogX)];

% Compute the geometric Brownian motion
X = exp(logX);

% Compute the expected path of the GBM
EX = exp(mu*t);

% Plot the expected, mean and sample paths of the GBM
figure(3)
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
ylim([0,2.5])
title('Paths of a geometric Brownian motion dX = \muXdt + \sigmaXdW')
print('-dpdf','gbppaths.pdf')

% Plot the probability density function of the GBM at different times
figure(4)

[h,x] = hist(X(40,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,1)
bar(x,f)
ylabel('f_X(x,0.2)')
xlim([0,2.5])
ylim([0,3])
title('Probability density function of a geometric Brownian motion at different times')

[h,x] = hist(X(100,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,2)
bar(x,f)
xlim([0,2.5])
ylim([0,3])
ylabel('f_X(x,0.5)')

[h,x] = hist(X(end,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,3)
bar(x,f)
xlim([0,2.5])
ylim([0,3])
xlabel('x')
ylabel('f_X(x,1)')

print('-dpdf','gbpdensities.pdf')