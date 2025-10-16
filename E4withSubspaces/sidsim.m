function [vPhi, vE, vQ] = sidsim(theta,din,T, i, n, s, ext)
%
%   [vPhi, vE, vQ] = sidsim(theta,din,T, i, n, s, ext)
%

if nargin < 7, ext = 0; end
if nargin < 6,   s = 1; end

N = 200;
m = din(1,2);

vPhi = zeros(n*m,N);
vE   = zeros(n*m,N);
vQ   = zeros(m*m,N);

for j=1:N
%
    fprintf(1,'%i',j);
    y = simgarch(theta,din,T);
%    y = simmod(theta,din,T);
y = y.^2 - mean(y.^2);
    [Phi,H,E,Q] = sident(y, [], i, n, s);
    [Phi, H, ign, E] = arma(Phi, H, E);
    F = Phi(:,1:m);

    vPhi(:,j) = F(:);
    vE(:,j)   = E(:);
    vQ(:,j)   = Q(:);
%
end

vPhi = sort(vPhi');
vE = sort(vE');
vQ = sort(vQ');