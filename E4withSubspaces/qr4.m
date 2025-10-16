function [HH, QZ, QT, r, ldetR12, M] = qr4(H)
% 7/6/07
% Copyright (c) Casals, Jerez and Sotoca, 2007

global E4OPTION

tol = E4OPTION(15);
m = size(H,1);
n = size(H,2);
r0 = min(m,n);

[Q,R,E] = qr(H');

if r0 > 1
   r = sum(abs(diag(R)) > tol);
else
   r = (abs(R(1,1)) > tol)*1;
end

QZ = Q(:,1:r)';
QT = Q(:,r+1:n)';
HH = E*R(1:r,1:m)';
ldetR12 = sum(log(diag(R(1:r,1:r).^2)));
M = E*[ones(r,1);zeros(m-r,1)] == 1;