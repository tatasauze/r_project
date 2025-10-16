function f = preesti3(theta, din, R, ix, ij, userf)
%
% 7/3/97
% Copyright (c) Jaime Terceiro, 1997

[Phi, Gam, E, H, D, C, Q, S, Rm] = thd2ee(theta, din, userf);

n = size(Phi,1);
m = size(H,1);

i = ij(1);
j = ij(2);
pond = ij(3);

O = [H; zeros((i-1)*m,n)];

Phi2 = Phi;

for j=1:i-1
    O(j*m+1:(j+1)*m,:) = H*Phi2;
    Phi2 = Phi2*Phi;
end

CRCt = C*Rm*C'; ESCt = E*S*C'; EQEt = E*Q*E';
P0 = lyapunov(Phi,EQEt);

for j=1:i+1
    B1  = H*P0*H' + CRCt;
    U = cholp(B1);
    K1  = ((Phi*P0*H' + ESCt)/U)/U';
    Phib= Phi - K1*H;
    P0  = Phib*P0*Phib' + EQEt + K1*CRCt*K1' - K1*ESCt' - ESCt*K1';
end

j=i;

iO = pinv(O);
O1 = O(1:(i-1)*m,:);

His = [[zeros(m,m);O(1:(j-2)*m,:)*E] zeros((j-1)*m,(j-2)*m)];
for k=2:j-1
    His(k*m+1:(j-1)*m,(k-1)*m+1:k*m) = His(m+1:(j-k)*m,1:m);
end
%kron(eye(j-1),CRCt) + His*kron(eye(j-1),Q)*His' + His*kron(eye(j-1),S*C') + (His*kron(eye(j-1),S*C')'
%O1*P0*O1'

V = cholp(O1*P0*O1' + kron(eye(j-1),CRCt) + His*kron(eye(j-1),Q)*His' + His*kron(eye(j-1),S*C') + (His*kron(eye(j-1),S*C'))');
%P0
%V'*V
%R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))'

X = [iO*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)) zeros(n,m)];

%e = V' \ [R(ix(3,1):ix(3,2),ix(1,1):ix(2,2))-O1*((Phi-E*H)*X+E*R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))), R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
e = [R(ix(3,1):ix(3,2),ix(1,1):ix(2,2))-O1*((Phi-E*H)*X+E*R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))), R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
f = log(det((e*e')));

%f = trace(e*e') + 2*sum(log(diag(V)));