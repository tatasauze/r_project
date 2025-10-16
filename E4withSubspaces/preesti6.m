function f = preesti6(theta, din, R, ij, userf)
%
% 7/3/97
% Copyright (c) Jaime Terceiro, 1997

[Phi, Gam, E, H, D, C, Q, S, Rm] = thd2ee(theta, din, userf);

n = size(Phi,1);
m = size(H,1);
r = max([size(Gam,2), size(D,2)]);

O = [H; zeros((ij-1)*m,n)];

Phi2 = Phi;

for k=1:ij-1
    O(k*m+1:(k+1)*m,:) = H*Phi2;
    Phi2 = Phi2*Phi;
end

CRCt = C*Rm*C'; SCt = S*C'; EQEt = E*Q*E';

P0 = lyapunov(Phi,EQEt);
His = [[zeros(m,m);O(1:(ij-1)*m,:)*E] zeros(ij*m,(ij-1)*m)];

for k=2:ij
    His(k*m+1:ij*m,(k-1)*m+1:k*m) = His(m+1:(ij-k+1)*m,1:m);
end

Iij = eye(ij);
P1 = His*kron(eye(ij),SCt);
V = cholp(O*P0*O' + kron(Iij,CRCt) + His*kron(Iij,Q)*His' + P1 + P1');

if r
   Op = eye(ij*m)- O*pinv(O);
   Hid = [[D;O(1:(ij-1)*m,:)*Gam] zeros(ij*m,(ij-1)*r)];
   for k=2:ij
       Hid((k-1)*m+1:ij*m,(k-1)*r+1:k*r) = Hid(1:(ij-k+1)*m,1:r);
   end
   e = V' \ [Op*(R(ij*r+1:ij*(m+r),1:ij*r)-Hid*R(1:ij*r,1:ij*r)) R(ij*r+1:ij*(m+r),ij*r+1:ij*(m+r))];
else
   e = V' \ R;
end

f = trace(e*e') + 2*sum(log(diag(V)));