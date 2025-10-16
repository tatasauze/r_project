function [vPhi, vsPhi, vE, vsE, vQ, vD, vsD, vD2] = simint(theta,din,u,T, N, i, ki, s, ext)
%
%  Preparado para simular únicamente modelos puramente estocásticos.
%  Simula para el rango i con método exacto y aprox., no ponderación y
%  ponderación var y O(i-1) estimada y pseudoinversa 
%  T : tamaño de la muestra
%  N : Nº de realizaciones
%  i : vector con los posibles valores de i
%  ki: índices de Kronecker (necesario especificar para univariantes)
%  s : periodo estacional (no sirve para nada)
% ext: matrices block-Hankel extendidas
%
%
% [vPhi, vE, vQ, vK] = simtot(theta,din,T, N, i, n, s, ext)

if nargin < 9, ext = 0; end
if nargin < 8,   s = 1; end
i = i(:);

n = sum(ki);

[Phi,Gam,E,H,D] = thd2ee(theta,din);
m = size(H,1);
r = size(D,2);

si = size(i,1);

vPhi = zeros(N,n*m);
vsPhi = zeros(N,n*m);
vE   = zeros(N,n*m);
vsE   = zeros(N,n*m);
vD   = zeros(N,m*r);
vD2   = zeros(N,m*r);
vsD   = zeros(N,m*r);
vQ   = zeros(N,m*m);

sidopt('verbose','no');

for j=1:N
%
    y = simmod(theta,din,T,u);
    vD2(j, :)   = (u \ y)';
    fmax = inf;
    for k=1:size(i,1)
        [Phi,sPhi,H,sH,E,sE,Q,D,sD,innov] = sidint(y, u, i(k), ki, s);
        [thet, di] = ee2thd(Phi, zeros(n,r), E+Phi(:,1), [H zeros(m,n-size(H,2))], D, [], Q);
        f = pem(thet,di,[y u]);
        if f < fmax
           fmax = f;
           Phix = Phi;
           Hx   = H;
           Ex   = E;
           Qx   = Q;
           Dx   = D;
           sPhix = sPhi;
           sHx   = sH;
           sEx   = sE;
           sDx   = sD;
        end
    end

    F = Phix(:,1:m);
    vPhi(j, :) = F(:)';
    vD(j, :)   = Dx(:)';
    vE(j, :)   = Ex(:)';
    vQ(j, :)   = Qx(:)';
    vsPhi(j, :) = sPhix(:)';
    vsD(j, :)   = sDx(:)';
    vsE(j, :)   = sEx(:)';
end
sidopt('verbose','si');
