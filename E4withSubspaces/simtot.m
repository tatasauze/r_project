function [vPhi, vE, vQ, vK] = simtot(theta,din,T, N, i, ki, s, ext)
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

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
i = i(:);

n = sum(ki);

[Phi,ign,E,H] = thd2ee(theta,din);
m = size(H,1);

si = size(i,1);

vPhi = zeros(N,n*m*si*8);
vE   = zeros(N,n*m*si*8);
vQ   = zeros(N,m*m*si*8);
vK   = zeros(N, si*8);

met = ['exa';'apr'];
pon = ['no';'va'];
obs = ['r';'m'];
sidopt('verbose','no');

for j=1:N
%
j
    y = simmod(theta,din,T);
    fmax = inf;
    ptr = 0;
    ptr2 = 0;
    ptr3 = 1;
    for k=1:size(i,1)
        for k1=1:2
            sidopt('met',met(k1,:));
        for k2=1:2
            sidopt('obs',obs(k2,:));
            for k3=1:2
                sidopt('pon',pon(k3,:));
                [Phi,H,E,Q,innov] = sident(y, [], i(k), n, s);
                Q = real(Q);
                [thet, di] = ee2thd(Phi, [], E, H, [], [], Q);
                f = fvmod(thet,di,y);
                [Phi, H, ign, E] = echelon(ki,Phi, H, E);
                F = Phi(:,1:m);

                vPhi(j, ptr+1:ptr+n*m) = F(:)';
                vE(j, ptr+1:ptr+n*m)   = E(:)';
                vQ(j, ptr2+1:ptr2+m*m)   = Q(:)';
                vK(j, ptr3) = f;
                ptr = ptr +n*m;
                ptr2 = ptr2 + m*m;
                ptr3 = ptr3 + 1;
            end
        end
        end
    end
end
sidopt('verbose','si');
