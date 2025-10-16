function [vPhi, vE, vQ] = simula(theta,din,T, N, i, n, s, ext)
%
% function [vPhi, vE, vQ] = simula(theta,din,T, N, i, n, s, ext)

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
i = i(:);

[Phi,ign,E,H] = thd2ee(theta,din);
m = size(H,1);

vPhi = zeros(n*m,N);
vE   = zeros(n*m,N);
vQ   = zeros(m*m,N);
vK   = zeros(N,1);

for j=1:N
%
    y = simmod(theta,din,T);
    fmax = inf;
    for k=i(1):i(size(i,1))
        [Phi,H,E,Q,innov] = sident(y, [], k, n, s);
        [thet, di] = ee2thd(Phi, [], E, H, [], [], Q);
        f = fvarmax(thet,di,y);
        if f < fmax
           kmax = k;
           Phimax = Phi;
           Hmax = H;
           Emax = E;
           Qmax = Q;
           fmax = f;
        end
    end
    [Phi, H, ign, E] = arma(Phimax, Hmax, Emax);
    F = Phi(:,1:m);

    vPhi(:,j) = F(:);
    vE(:,j)   = E(:);
    vQ(:,j)   = Q(:);
    vK(j)     = kmax;
%
end

%vPhi = sort(vPhi');
%vE = sort(vE');
%vQ = sort(vQ');