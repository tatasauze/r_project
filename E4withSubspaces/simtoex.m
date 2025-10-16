function [vPhi, vE, vQ, vD, vGam, vK] = simtot(theta,din,thetau,dinu,T, N, i, ki, s, ext)
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

if nargin < 10, ext = 0; end
if nargin < 9,   s = 1; end
i = i(:);

n = sum(ki);

[Phi,Gam,E,H,D] = thd2ee(theta,din);
m = size(H,1);
r = size(Gam,2);

si = size(i,1);

vPhi = zeros(N,n*m*si*8);
vGam = zeros(N,n*r*si*8);
vE   = zeros(N,n*m*si*8);
vD   = zeros(N,m*r*si*8);
vQ   = zeros(N,m*m*si*8);
vK   = zeros(N, si*8);

met = ['exa';'apr'];
pon = ['no';'va'];
obs = ['r';'m'];
sidopt('verbose','no');

for j=1:N
%
j
    u = simmod(thetau,dinu,T);
    y = simmod(theta,din,T,u);
    fmax = inf;
    ptr = 0;
    ptr2 = 0;
    ptr4 = 0;
    ptr5 = 0;
    ptr3 = 1;
    for k=1:size(i,1)
        for k1=1:2
            sidopt('met',met(k1,:));
        for k2=1:2
            sidopt('obs',obs(k2,:));
            for k3=1:2
                sidopt('pon',pon(k3,:));
                [Phi,H,E,Q,Gam,D,innov] = sident(y, u, i(k), n, s);
                Q = real(Q);
                [thet, di] = ee2thd(Phi, Gam, E, H, D, [], Q);
                f =  log(det(innov*innov'/size(innov,2)));

%                f = fvmod(thet,di,[y u]);
                [Phi, H, ign, E, G] = echelon(ki,Phi, H, E, Gam, D);
                F = Phi(:,1:m);

                vPhi(j, ptr+1:ptr+n*m) = F(:)';
                vGam(j, ptr4+1:ptr4+n*r) = G(:)';
                vD(j, ptr5+1:ptr5+m*r) = D(:)';
                vE(j, ptr+1:ptr+n*m)   = E(:)';
                vQ(j, ptr2+1:ptr2+m*m)   = Q(:)';
                vK(j, ptr3) = f;
                ptr = ptr +n*m;
                ptr2 = ptr2 + m*m;
                ptr3 = ptr3 + 1;
                ptr4 = ptr4 + n*r;
                ptr5 = ptr5 + m*r;
            end
        end
        end
    end
end
sidopt('verbose','si');
