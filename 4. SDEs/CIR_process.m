%% Simulate a Cox-Ingersoll-Ross process
% dX = alpha*(mu-X)*dt + sigma*sqrt(X)*dW

% Define parameters and time grid
npaths = 20000; % number of paths
T = 1; % time horizon
nsteps = 200; % number of time steps
dt = T/nsteps; % time step
t = 0:dt:T; % observation times
alpha = 5; mu = 0.07; sigma = 0.265; % model parameters
X0 = 0.03; % initial value

%% Allocate and initialise all paths
X = [X0*ones(1,npaths);zeros(nsteps,npaths)];

% Sample standard Gaussian random numbers
N = randn(nsteps,npaths);

% Degrees of freedom of the non-central chi square distribution
%d = 4*alpha*mu/sigma^2;

% Compute and accumulate the increments
a = sigma^2/alpha*(exp(-alpha*dt)-exp(-2*alpha*dt)); % Euler with analytic moments
b = mu*sigma^2/(2*alpha)*(1-exp(-alpha*dt))^2; % Euler with analytic moments
%k = sigma^2*(1-exp(-alpha*dt))/(4*alpha); % exact method
for i = 1:nsteps
    %X(i+1,:) = X(i,:) + alpha*(mu-X(i,:))*dt + sigma*sqrt(X(i,:)*dt).*N(i,:); % plain Euler
    X(i+1,:) = mu+(X(i,:)-mu)*exp(-alpha*dt) + sqrt(a*X(i,:)+b).*N(i,:); % Euler with a.m.
    %lambda = 4*alpha*X(i,:)/(sigma^2*(exp(alpha*dt)-1)); % exact method
    %X(i+1,:) = icdf('ncx2',rand(1,npaths),d,lambda)*k; % exact method
end

% Compute the expected path
EX = mu + (X0-mu)*exp(-alpha*t);

%% Plot the expected, mean and sample paths
figure(1)
plot(t,EX,'k',t,mean(X,2),':k',t,X(:,1:1000:end),t,EX,'k',t,mean(X,2),':k')
legend('Expected path','Mean path')
xlabel('t')
ylabel('X')
sdevinfty = sigma*sqrt(mu/(2*alpha));
ylim([0 mu+4*sdevinfty])
title('Paths of a Cox-Ingersoll-Ross process dX = \alpha(\mu-X)dt + \sigmaX^{1/2}dW')
print('-dpdf','cirppaths.pdf')

%% Compute and plot the probability density function at different times
t = [0.05 0.1 0.2 0.4 1];
x = linspace(0,mu+4*sdevinfty,200);
k = sigma^2*(1-exp(-alpha*t))/(4*alpha);
d = 4*alpha*mu/sigma^2;
lambda = 4*alpha*X0./(sigma^2*(exp(alpha*t)-1)); % non-centrality parameter
f = zeros(length(x),length(t));
for i = 1:length(t)
    f(:,i) = pdf('ncx2',x/k(i),d,lambda(i))/k(i);
end

figure(3)
plot(x,f)
xlabel('x')
ylabel('f_X(x,t)')
legend('t = 0.05','t = 0.10','t = 0.20','t = 0.40','t = 1.00')
title('Probability density function of a Cox-Ingersoll-Ross process at different times')
print('-dpdf','cirpdensities.pdf')