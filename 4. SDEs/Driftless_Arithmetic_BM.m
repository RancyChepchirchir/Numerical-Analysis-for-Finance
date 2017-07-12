% Einstein relation for a driftless arithmetic Brownian motion dX(t) = sigma*dW(t)

% Define parameters and time grid
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
sigma = 1; % model parameters

% Compute the increments
dX = sigma*sqrt(dt)*randn(nsteps,npaths);

% Accumulate the increments
X = [zeros(1,npaths); cumsum(dX)];

% Compute the expected path
EX = zeros(size(t));

% Plot the expected, mean and sample path
figure(1)
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
ylim([-2,2]);
title('Paths of a driftless arithmetic Brownian motion dX(t) = \sigmadW(t)')

% Plot the probability density function at different times
figure(2)

[h,x] = hist(X(40,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,1)
bar(x,f)
ylabel('f_X(x,0.2)')
xlim([-3,3])
ylim([0,1])
title('Probability density function of a driftless arithmetic Brownian motion at different times')

[h,x] = hist(X(100,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,2)
bar(x,f)
xlim([-3,3])
ylim([0,1])
ylabel('f_X(x,0.5)')

[h,x] = hist(X(end,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,3)
bar(x,f)
xlim([-3,3])
ylim([0,1])
xlabel('x')
ylabel('f_X(x,1)')

% Plot the mean square displacement
figure(3)
msd = mean(X.^2,2);
sigsq = mean(msd(2:end)'./t(2:end));
sig = sqrt(sigsq)
plot(t,msd,'r',t,sigsq*t,'g')
legend('Sampled','\sigma^2t')
xlabel('t')
ylabel('E[(X-X_0)^2]')
title('Einstein relation for dX(t) = \sigmadW(t)')
%print('-dpdf','msd.pdf')