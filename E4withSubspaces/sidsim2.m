function [vPhi, vE, vQ, vD, vG] = sidsim(theta,din,T, i, n, s, ext)
%
%   [vPhi, vE, vQ] = sidsim(theta,din,T, i, n, s, ext)
%

if nargin < 7, ext = 0; end
if nargin < 6,   s = 1; end

N = 500;
m = din(1,2);
r = din(1,3);

vPhi = zeros(n*m*2,N);
vE   = zeros(n*m*2,N);
vQ   = zeros(m*m*2,N);
vD   = zeros(m*r,N);
vG   = zeros(n*r,N);

for j=1:N
%
    fprintf(1,'%i',j);
%    y = simgarch(theta,din,T);
    u = zeros(T,1);
    u(T/2) = 1;
    y = simmod(theta,din,T,u);
%y = y.^2 - mean(y.^2);
    [Phi,H,E,Q,Gam,D] = sident(y, u, i, n, s);
    [Phi, H, ign, E, G] = arma(Phi, H, E,Gam,D);
    F = Phi(:,1:m);

    vPhi(1:n,j) = F(:);
    vE(1:n,j)   = E(:);
    vQ(1:m,j)   = Q(:);
    vG(:,j)   = G(:);
    vD(:,j)   = D(:);

    [Phi,H,E,Q] = sident(y, [], i, n, s);
    [Phi, H, ign, E] = arma(Phi, H, E);
 
    F = Phi(:,1:m);

    vPhi(n+1:2*n,j) = F(:);
    vE(n+1:2*n,j)   = E(:);
    vQ(m+1:2*n,j)   = Q(:);

%
end

vPhi = sort(vPhi');
vE = sort(vE');
vQ = sort(vQ');
vG = sort(vG');
vD = sort(vD');