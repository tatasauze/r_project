function [vPhi, vE, vQ1, vQ2, vR] = simerv(theta,din,T, N, i, ki, s, ext)
%
%  Simula modelos estacionales con matrices block-Hankel extendidas y no extendidas
%  Atención, fuera de la función es preciso poner las opciones de estimación por 
%  defecto
%
% [vPhi, vE, vQ, vK] = simtot(theta,din,T, N, i, n, s, ext)

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
i = i(:);

n = sum(ki);

[Phi,ign,E,H] = thd2ss(theta,din);
m = size(H,1);

vPhi = zeros(N,n*m);
vE   = zeros(N,n*m);
vQ1  = zeros(N,m*m);
vQ2  = zeros(N,m*m);
vR   = zeros(N,m*m);

sidopt('verbose','no');

for j=1:N
%
j
    y = simmod(theta,din,T);
    fmax = inf;
    for k=1:size(i,1)
        [Phi,H,E,Q,innov,G,L0] = sident(y, [], i(k), n, s);
        [thet, di] = ss2thd(Phi, [], E, H, [], [], Q);
        f = pem(thet,di,y);
%        if trace(innov*innov'/size(innov,2)) < fmax
%           fmax =  trace(innov*innov'/size(innov,2));
        if f < fmax
           fmax = f;
           imax = i(k);
           Phix = Phi;
           Hx   = H;
           Ex   = E;
           Qx   = Q;
           Gx   = G;
           L0x  = L0;
        end
    end
    [Phi, H, Tr, E] = echelon(ki,Phix, Hx, Ex);
    G = inv(Tr)*Gx;
    F = Phi(:,1:m);
    vPhi(j,:) = F(:)';
    vE(j,:)   = E(:)';
    [P, Q2, R] = siderv(Phi, H, F, G, L0x);
    vQ1(j, :)   = Qx(:)';
    vQ2(j, :)   = Q2(:)';
    vR(j, :)   = R(:)';
end
sidopt('verbose','si');
