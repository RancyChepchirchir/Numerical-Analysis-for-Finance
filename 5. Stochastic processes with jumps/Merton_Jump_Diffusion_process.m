% Simulate a Merton jump-diffusion process
% X(t) = mu*t + sigma*W(t) + sum_{i=1}^{N(t)} Z_i 

% Define parameters and time grid
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
mu = 0.2; sigma = 0.3; % model parameters of the diffusion part (ABM)
muZ = 0.1; sigmaZ = 0.15; lambda = 0.5; % model parameters of the jump part (NCPP)

% Compute the increments of the ABM
dW = mu*dt + sigma*sqrt(dt)*randn(nsteps,npaths);

% Compute the increments of the NCPP
dN = poissrnd(lambda*dt,[nsteps,npaths]);
dJ = muZ*dN + sigmaZ*sqrt(dN).*randn(nsteps,npaths);

% Sum the increments of the ABM and the NCPP
dX = dW + dJ;

% Accumulate the increments
X = [zeros(1,npaths); cumsum(dX)];

% Compute the expected path
EX = (mu+lambda*muZ)*t;

% Plot the expected, mean and sample path
figure(3)
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
ylim([-0.8,1.2]);
title('Paths of a Merton jump-diffusion process X = \mut + \sigmaW(t) + \Sigma_{i=1}^{N(t)} Z_i')
print('-dpdf','mjdpaths.pdf')

% Plot the probability density function at different times
figure(4)

[h,x] = hist(X(40,:),100);
f = h/(sum(h)*(x(2)-x(1)));
subplot(3,1,1)
bar(x,f)
ylabel('f_X(x,0.2)')
xlim([-0.8,1.2])
ylim([0,3])
title('Probability density function of a Merton jump-diffusion process at different times')

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

print('-dpdf','mjddensities.pdf')