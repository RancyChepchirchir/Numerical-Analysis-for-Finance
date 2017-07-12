%% Compute the Black-Scholes price of a European option
% The model for the underlying is geometric Brownian motion
% dS = mu*S*dt + sigma*S*dW

% Model parameters
T = 1; % maturity
S0 = 1; % spot price
K = 1.1; % strike price
r = 0.05; % risk-free interest rate
q = 0.02; % dividend rate
sigma = 0.4; % volatility
mu = r-q-0.5*sigma^2; % risk-neutral measure

% Monte Carlo parameters
nblocks = 2000; % number of blocks
nsample = 100000; % number of samples per block

% Fourier parameters
xwidth = 6; % width of the support in real space
ngrid = 2^8; % number of grid points
alpha = -10; % damping factor for call

%% Analytical solution
tic
d1 = (log(S0/K)+(r-q+0.5*sigma^2)*T)/(sigma*sqrt(T));
% d2 = d1 - sigma*sqrt(T);
d2 = (log(S0/K)+(r-q-0.5*sigma^2)*T)/(sigma*sqrt(T));
Vca = S0*exp(-q*T)*cdf('norm',d1,0,1) - K*exp(-r*T)*cdf('norm',d2,0,1);
% Put-call parity: V_call + K*exp(-r*T) = V_put - S_0
Vpa = K*exp(-r*T)*cdf('norm',-d2,0,1) - S0*exp(-q*T)*cdf('norm',-d1,0,1);
cputime_a = toc;

% Analytical solution provided by Matlab's Financial Toolbox
tic
[VcaM,VpaM] = blsprice(S0,K,r,T,sigma,q);
cputime_aM = toc;

% Print the results
fprintf('%17s%14s%14s%14s\n','','call','put','CPU_time/s')
fprintf('%17s%14.10f%14.10f%14.10f\n','Analytical',Vca,Vpa,cputime_a)
fprintf('%17s%14.10f%14.10f%14.10f\n','Analytical Matlab',VcaM,VpaM,cputime_aM)

% Plot the analytical solution
[Sa,t] = meshgrid(0:.05:2,0:0.025:T);
d1 = (log(Sa/K)+(r-q+0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
d2 = (log(Sa/K)+(r-q-0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));

close all
figure(1)
Vc = Sa.*exp(-q*(T-t)).*cdf('norm',d1,0,1) - K*exp(-r*(T-t)).*cdf('norm',d2,0,1);
Vc(end,:) = max(Sa(end,:)-K,0);
mesh(Sa,t,Vc)
xlabel('S')
ylabel('t')
zlabel('V')
title('Call')
view(-30,24)
print('-dpng','bsc.png')

figure(2)
Vp = K*exp(-r*(T-t)).*cdf('norm',-d2,0,1) - Sa.*exp(-q*(T-t)).*cdf('norm',-d1,0,1);
Vp(end,:) = max(K-Sa(end,:),0);
mesh(Sa,t,Vp)
xlabel('S')
ylabel('t')
zlabel('V')
title('Put')
view(30,24)
print('-dpng','bsp.png')

% Plot the analytical solution as a function of the log price
k = log(K/S0);
[xa,t] = meshgrid(-1:.05:1,0:0.025:T);
d1 = (xa-k+(r-q+0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));
d2 = (xa-k+(r-q-0.5*sigma^2)*(T-t))./(sigma*sqrt(T-t));

figure(3)
Vc = S0*(exp(xa-q*(T-t)).*cdf('norm',d1,0,1) - exp(k-r*(T-t)).*cdf('norm',d2,0,1));
Vc(end,:) = S0*max(exp(xa(end,:))-exp(k),0);
mesh(xa,t,Vc)
xlabel('x')
ylabel('t')
zlabel('V')
title('Call')
view(-30,24)
print('-dpng','bscx.png')

figure(4)
Vp = S0*(exp(k-r*(T-t)).*cdf('norm',-d2,0,1) - exp(xa-q*(T-t)).*cdf('norm',-d1,0,1));
Vp(end,:) = S0*max(exp(k)-exp(xa(end,:)),0);
mesh(xa,t,Vp)
xlabel('x')
ylabel('t')
zlabel('V')
title('Put')
view(30,24)
print('-dpng','bspx.png')

%% Fourier method

% Grids in real and Fourier space
tic
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid;
x = dx*(-N:N-1);
dxi = 2*pi/xwidth; % Nyquist relation
xi = dxi*(-N:N-1);

% Characteristic function at time T
xia = xi+1i*alpha; % call
psi = 1i*mu*xia-0.5*(sigma*xia).^2; % characteristic exponent
Psic = exp(psi*T); % characteristic function
xia = xi-1i*alpha; % put
psi = 1i*mu*xia-0.5*(sigma*xia).^2; % characteristic exponent
Psip = exp(psi*T); % characteristic function

% These functions provide the characteristic functions of 8 Levy processes
% param = parameters(1,T,T,r,q); % sigma must be set inside parameters.m
% [x,fc,xi,Psic] = kernel(ngrid,-b,b,param,alpha,0,1); % call
% [x,fp,xi,Psip] = kernel(ngrid,-b,b,param,-alpha,0,1); % put

% Fourier transform of the payoff
U = S0*exp(b);
L = S0*exp(-b);
[~,gc,Gc] = payoff(x,xi,alpha,K,L,U,S0,1); % call
[S,gp,Gp] = payoff(x,xi,-alpha,K,L,U,S0,0); % put

% Discounted expected payoff computed with the Plancherel theorem
c = exp(-r*T).*real(fftshift(fft(ifftshift(Gc.*conj(Psic)))))/xwidth; % call
VcF = interp1(S,c,S0,'spline');
p = exp(-r*T).*real(fftshift(fft(ifftshift(Gp.*conj(Psip)))))/xwidth; % put
VpF = interp1(S,p,S0,'spline');
cputime_F = toc;
fprintf('%17s%14.10f%14.10f%14.10f\n','Fourier',VcF,VpF,cputime_F)

% Figures
% figures_ft(S,x,xi,Psic,gc,Gc) % call
% figures_ft(S,x,xi,Psip,gp,Gp) % put

%% Monte Carlo

tic;
VcMCb = zeros(nblocks,1);
VpMCb = zeros(nblocks,1);
for i = 1:nblocks

    % Arithmetic Brownian motion X = log S at time T
    X = mu*T + sigma*sqrt(T)*randn(1,nsample);

    % Geometric Brownian motion at time T
    S = exp(X);

    % Discounted expected payoff
    VcMCb(i) = exp(-r*T)*mean(max(S-K,0));
    VpMCb(i) = exp(-r*T)*mean(max(K-S,0));

end
VcMC = mean(VcMCb);
VpMC = mean(VpMCb);
scMC = sqrt(var(VcMCb)/nblocks);
spMC = sqrt(var(VpMCb)/nblocks);
cputime_MC = toc;
fprintf('%17s%14.10f%14.10f%14.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)
fprintf('%17s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)
