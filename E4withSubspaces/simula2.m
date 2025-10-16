function [vPhi, vE, vQ, vK, vsPhi, vsE] = simula2(theta,din,T, N, i, n, s, ext)
%
% function [vPhi, vE, vQ, vK] = simula(theta,din,T, N, i, n, s, ext)

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
i = i(:);

[Phi,ign,E,H] = thd2ee(theta,din);
m = size(H,1);

vPhi = zeros(n*m,N);
vE   = zeros(n*m,N);
vQ   = zeros(m*m,N);
vK   = zeros(N,1);
vsPhi = zeros(n*m,N);
vsE   = zeros(n*m,N);

for j=1:N
%
fprintf(1,'%1d ',j)
    y = simmod(theta,din,T);
    fmax = inf;
    for k=i(1):i(size(i,1))
        [Phi,sPhi,H,sH,E,sE,Q] = sideche3(y, [], k, n, s);
        Q = abs(Q);
        [thet, di] = ee2thd(Phi, [], E, H, [], [], Q);
        f = fvarmax(thet,di,y);
        if f < fmax
           kmax = k;
           Phimax = Phi;
           sPhimax = sPhi;
           Emax = E;
           sEmax = sE;
           Qmax = Q;
           fmax = f;
        end
    end
    F = Phimax(:,1:m);
    sF = sPhimax(:,1:m);
    vPhi(:,j) = F(:);
    vE(:,j)   = Emax(:);
    vsPhi(:,j) = sF(:);
    vsE(:,j)   = sEmax(:);
    vQ(:,j)   = Qmax(:);
    vK(j)     = kmax;
%
end

%vPhi = sort(vPhi');
%vE = sort(vE');
%vQ = sort(vQ');