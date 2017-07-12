% Plot the non-central chi-square distribution

clear all;

% Assign parameters
n = 5; % degrees of freedom
noncentralpar = 2; % non-centrality parameter

% Assign x range
x = linspace(0,20,100);

% Compute PDF and CFF
pdfncchi2 = pdf('ncx2',x,n,noncentralpar);
cdfncchi2 = cdf('ncx2',x,n,noncentralpar);

% Make the plot
plot(x,pdfncchi2,'r',x,cdfncchi2,'b')
xlabel('x')
title('Non central Chi-Square pdf and cdf when n=5 and d=2')
legend('Non Central Chi-Square pdf','Non Central Chi-Square cdf')
print('-dpdf','ncchi2.pdf')