% Compute the scale, the payoff and its Fourier transform
function [S,g,G] = payoff(x,xi,alpha,K,L,U,C,call)

% Scale
S = C*exp(x);

% Payoff 
if call == 1 % see e.g. Green, Fusai, Abrahams 2010, Eq. (3.24)
    g = exp(alpha*x).*max(S-K,0).*(S>=L).*(S<=U);
else % put
    g = exp(alpha*x).*max(K-S,0).*(S>=L).*(S<=U);
end

% Analytical Fourier transform of the payoff
l = log(L/C);
k = log(K/C);
u = log(U/C);

% Integration bounds
if call == 1 % call
    a = max(l,k);
    b = u;
else % put
    a = min(k,u);
    b = l;
end

% Green, Fusai, Abrahams 2010 Eq. (3.26) with extension to put option
xi2 = alpha+1i*xi;
G = C*((exp(b*(1+xi2))-exp(a*(1+xi2)))./(1+xi2) ...
    - (exp(k+b*xi2)-exp(k+a*xi2))./xi2);

% Eliminable discontinuities for xi = 0, otherwise 0/0 = NaN
if (alpha == 0)
    G(floor(end/2)+1) = C*(exp(b)-exp(a)-exp(k)*(b-a));
elseif (alpha == -1)
    G(floor(end/2)+1) = C*(b-a+exp(k-b)-exp(k-a));
end
