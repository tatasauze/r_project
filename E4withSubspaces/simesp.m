function [vPhi, vE, vQ, vGam, vD, vK, vi, vPK] = simesp(theta,din,thetau,dinu, T, N, i, ki, option, optval)
%
%  Especificación automática de modelos
%  ATENCION, POR EL MOMENTO NO IDENTIFICA CON INDICES DE KRONECKER
%
%  [vPhi, vE, vQ, vGam, vD, vK, vi] = simesp(theta,din,thetau,dinu, T, N, i, ki, {option, optval})

i = i(:);
ix = size(i,1);

if nargin < 10, option = []; optval = 1; end
sopt = size(optval,1);

n = sum(ki);

[Phi, Gam, E, H, D, C, Q] = thd2ee(theta, din);
m = size(H,1);
r = max(size(Gam,2),size(D,2));

kk = n/m;

vPhi = zeros(N,n*m*sopt);
vE   = zeros(N,n*m*sopt);
vQ   = zeros(N,m*m*sopt);
vK   = zeros(N,1);
vi   = zeros(N,sopt);
vGam = zeros(N,n*r*sopt);
vD   = zeros(N,m*r*sopt);
vPK  = zeros(N,2*m*n);

sidopt('verbose','no');

for j=1:N
%
    if r
       u = simmod(thetau, dinu,T);
       y = simmod(theta,din,T,u);
    else
       u = [];
       y = simmod(theta,din,T);
    end

    t = sidang(y, u, i(ix));
    kix = find(diag(t) < .95);

    if isempty(kix), kix = i(ix)*m; else kix = kix(1)-1; end
    vK(j) = kix;

    for k1=1:sopt
      if sopt > 1, sidopt(option,optval(k1,:)); end
    
      fmax = inf;
      for k=1:size(i,1)
          if r
             [Phi,H,E,Q,Gam,D] = sident(y, u, i(k), n);
             [thet, di] = ee2thd(Phi, Gam, E, H, D, [], Q);
             f = pem(thet,di,[y u]);
          else
             [Phi,H,E,Q] = sident(y, [], i(k), n);
             [thet, di] = ee2thd(Phi, [], E, H, [], [], Q);
             f = pem(thet,di,y);
          end
          if f < fmax
             fmax = f;
             imax = i(k);
             Phix = Phi;
             Hx   = H;
             Ex   = E;
             Qx   = Q;
             if r
                Gamx   = Gam;
                Dx   = D;
             end
          end
      end
      if r
         [Phi, H, Tr, E, G] = echelon(ki,Phix, Hx, Ex, Gamx, Dx);
         vGam(j,(k1-1)*n*r+1:k1*n*r) = G(:)';
         vD(j,(k1-1)*m*r+1:k1*m*r) = Dx(:)';
      else
         [Phi, H, Tr, E] = echelon(ki,Phix, Hx, Ex);
      end
      F = Phi(:,1:m);
      vPhi(j,(k1-1)*n*m+1:k1*n*m) = F(:)';
      vE(j,(k1-1)*n*m+1:k1*n*m)   = E(:)';
      vQ(j,(k1-1)*m*m+1:k1*m*m)   = Qx(:)';
      vi(j,k1) = imax;

    end
    if nargin > 7
       vPK(j,:) = pk(y,kk,kk);
    end
end
sidopt('verbose','si');
