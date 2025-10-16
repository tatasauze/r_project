function [vPhi, vPhi2, vE, vE2, vD, vD2, vGam, vGam2] = simula(thetay,diny,thetau,dinu,T, N, i, n, s, ext)
%
% function [vPhi, vE, vQ] = simula(theta,din,T, N, i, n, s, ext)

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
i = i(:);

[Phi,Gam,E,H] = thd2ee(thetay,diny);
m = size(H,1);
r = size(Gam,2);

vPhi = zeros(n*m,N);
vE   = zeros(n*m,N);
vPhi2 = zeros(n*m,N);
vE2   = zeros(n*m,N);
vD    = zeros(m*r,N);
vD2   = zeros(m*r,N);
vGam  = zeros(n*r,N);
vGam2 = zeros(n*r,N);

for j=1:N
%
j
    u = simmod(thetau,dinu,T);
    y = simmod(thetay,diny,T,u);
    [TH A0] = n4sid2([y u],n,size(y,2),i,[1 1 0]);
    [A B C D K] = th2ss(TH);
    [Phi, H, Tar, E, Gam] = arma(A, C, K, B, D);
    [Phi2,H2,E2,Q,Gam2,D2] = sidentpr(y, u, i, n);
    [Phi2, H2, Tar, E2, Gam2] = arma(Phi2, H2, E2, Gam2, D2);

    F = Phi(:,1:m);
    vPhi(:,j) = F(:);
    vE(:,j)   = E(:);
    vGam(:,j) = Gam(:);
    vD(:,j)   = D(:);
    F = Phi2(:,1:m);
    vPhi2(:,j) = F(:);
    vE2(:,j)   = E2(:);
    vGam2(:,j) = Gam2(:);
    vD2(:,j)   = D2(:);

%
end

vPhi = vPhi';
vE = vE';
vGam = vGam';
vD = vD';
vPhi2 = vPhi2';
vE2 = vE2';
vGam2 = vGam2';
vD2 = vD2';
