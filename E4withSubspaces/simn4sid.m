function [vPhi, vPhi2, vE, vE2] = simula(theta,din,T, N, i, n, s, ext)
%
% function [vPhi, vE, vQ] = simula(theta,din,T, N, i, n, s, ext)

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
i = i(:);

[Phi,ign,E,H] = thd2ee(theta,din);
m = size(H,1);

vPhi = zeros(n*m,N);
vE   = zeros(n*m,N);
vPhi2 = zeros(n*m,N);
vE2   = zeros(n*m,N);

for j=1:N
%
j
    y = simmod(theta,din,T);
    [TH A0] = n4sid(y,n,size(y,2),i);
    [A B C D K] = th2ss(TH);
    [Phi, H, Tar, E] = arma(A, C, K);
    [Phi2,H2,E2,Q,innov] = sident(y, [], i, n);
    [Phi2, H2, ign, E2] = arma(Phi2, H2, E2);

    F = Phi(:,1:m);
    vPhi(:,j) = F(:);
    vE(:,j)   = E(:);
    F = Phi2(:,1:m);
    vPhi2(:,j) = F(:);
    vE2(:,j)   = E2(:);

%
end

vPhi = vPhi';
vE = vE';
vPhi2 = vPhi2';
vE2 = vE2';
