% Linear congruential random number generator
% Seydel, Course Notes, Chapter 2, pages 202-203

%seed = 12345;
seed = tic/10^9;
a = 1597; b = 51749; M = 244944;
nsample = 10000;

% Generator
U = zeros(1,nsample);
U(1) = seed;
for i = 2:nsample
    U(i) = mod((a*U(i-1)+b),M);  
end
U = U/M;

% Bin the random variables in a histogram and normalise it
nbins = 100;
[h,x] = hist(U,100);
h = h*nbins/nsample; % normalisation

% Plot of the probability density function
figure(1)
plot(x,h,x,ones(size(x)))
ylim([0 2])
xlabel('x')
ylabel('f')
legend('Sampled','Theory')
title('Uniform distribution')