function f=exgausspdf(mu,sig,tau,x)

% given parameter mu, sig and tau, returns density at x for the ex-Gaussian
% mu, sig, tau are scalars
% x is either a scalar, vector or matrix
% f has the same shape as x

arg1=(mu./tau)+((sig.*sig)./(2.*tau.*tau))-(x/tau);
arg2=((x-mu)-((sig.*sig)./tau))./sig;
f=(1./tau)*(exp(arg1).*pnf(arg2));

function p=pnf(x) 
% cumulative probability for the normalized Gaussian
% x is scalar, vector or matrix % f has the same shape as x 

a=find(x<0); b=find(x>=0); p=x; 
m_sqrt2=sqrt(2); 
p(b)=( (1+erf(x(b)./m_sqrt2)) ./ 2 ); 
p(a)=( (erfc(-x(a)./m_sqrt2)) ./ 2 );