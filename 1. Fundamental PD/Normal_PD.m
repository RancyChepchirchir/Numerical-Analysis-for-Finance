% Plot and sample the normal distribution

clear all
close all

% Define the parameters
mu = 0.2;
sigma = 0.1;
nsample = 10^6;

% Define the grid
%x = linspace(-.4,0.8,121);
x = -0.4:0.01:0.8;

% Compute the PDF and the CDF
%f = 1/(sqrt(2*pi)*sigma)*exp(-((x-mu)/sigma).^2/2);
f = pdf('norm',x,mu,sigma);
F = cdf('norm',x,mu,sigma);

% Open a plot window
figure(1)
plot(x,f,'r',x,F,'b');

% Define the limits of the x-axis
xlim([-0.4 0.8])

% Label the x axis
xlabel('x')

% Insert the legend
legend('PDF','CDF')

% Set the title
title('Normal distribution with \mu = 0.2 and \sigma = 0.1')

% Print the figure to a file
print('-dpdf','normal.pdf')

% Sample a normal distribution
%U = rand(nsample,1); % method 1
%X = mu + sigma*norminv(U,0,1); % method 1a
%X = norminv(U,mu,sigma); % method 1b 
X = mu + sigma*randn(nsample,1); % method 2

figure(2)
histfit(X,100)
legend('Sampled','Normal fit')
title('Normal distribution with \mu = 0.2 and \sigma = 0.1')

[h,xx] = hist(X,100);
h = h/(nsample*(xx(2)-xx(1)));

figure(3)
plot(xx,h,'b',x,f,'r')
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('Normal distribution with \mu = 0.2 and \sigma = 0.1')

figure(4)
bar(xx,h)
hold on
plot(x,f,'r')
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('Normal distribution with \mu = 0.2 and \sigma = 0.1')