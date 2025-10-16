function [theta,din,lab,std] = subest1seas(z,u,nn,s);

% SUBEST1SEAS - Computes a fast estimate of an univariate ARMA model
% [Seasonal or non-seasonal univariate series]
% [theta,din,lab,std] = subest1seas(z,u,[n ns],s)
%
% z         > time series
% u         > input series [],
% n         > regular subsystem order
% ns        > seasonal subsystem order
% s         > seasonality frequency
% theta     < vector of estimates
% din       < vector which stores a description of the model dynamics
% lab       < vector which contains the names for the parameters in theta
% std       < vector which contains the standard deviations of the parameters in theta

global E4OPTION

if nargin < 4, s = 1; end
n = nn(1); ns = nn(2); 
N = size(z,1);
r = size(u,2);
is = ns+3; 
if s == 4   
    i = n+1;
else
    i = max(round(log(N)),n+1); %% log(N)!!
end
if r
    [F,sF,TH,sTH,SIGMA,INNOV] = varmaech(z, [], is, ns, s);
    [f,sf,th,sth,sigma,g,sg] = varmaech(INNOV', u, i, n);
else
    [F,sF,TH,sTH,SIGMA,INNOV] = varmaech(z, [], is, ns, s);
    [f,sf,th,sth,sigma] = varmaech(INNOV', [], i, n);
    g = []; sg = []; 
end
[theta, din, lab] = arma2thd([f(2:n+1)], [F(2:ns+1)], [th(2:n+1)], [TH(2:ns+1)], sigma, s, g, r);
std = arma2thd([sf(2:n+1)], [sF(2:ns+1)], [sth(2:n+1)], [sTH(2:ns+1)], 0, s, sg, r);