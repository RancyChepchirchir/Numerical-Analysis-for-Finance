% Simulate a time-changed arithmetic Brownian motion: the variance Gamma process
% dX(t) = theta*dG(t) + sigma*dW(G(t))

% Define parameters and time grid
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
theta = 0.2; sigma = 0.3; kappa = 0.25; % model parameters

% Compute the increments of the Gamma process
dG = gamrnd(dt/kappa,kappa,[nsteps,npaths]);

% Compute the increments of the ABM on the Gamma random clock
dX = theta*dG + sigma*sqrt(dG).*randn(nsteps,npaths);

% Accumulate the increments
X = [zeros(1,npaths); cumsum(dX)];

% Compute the expected path
EX = theta*t;

% Plot the expected, mean and sample path
figure(3)
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
ylim([-0.8,1.2])
title('Paths of a variance Gamma process dX(t) = \mudG(t) + \sigmadW(G(t))')
print('-dpdf','vgpaths.pdf')

% Plot the probability density function at different times
figure(4)

[h,x] = hist(X(40,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,1)
bar(x,f)
ylabel('f_X(x,0.2)')
xlim([-0.8,1.2])
ylim([0,6])
title('Probability density function of a variance Gamma process at different times')

[h,x] = hist(X(100,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,2)
bar(x,f)
xlim([-0.8,1.2])
ylim([0,3])
ylabel('f_X(x,0.5)')

[h,x] = hist(X(end,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,3)
bar(x,f)
xlim([-0.8,1.2])
ylim([0,3])
xlabel('x')
ylabel('f_X(x,1)')

print('-dpdf','vgdensities.pdf')