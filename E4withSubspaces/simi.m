function [vPhi, vE, vQ, vGam, vD, vK] = simi(theta,din,thetau,dinu, T, N, i, ki, option, optval)
%
%  Especificación automática de modelos
%  ATENCION, POR EL MOMENTO NO IDENTIFICA CON INDICES DE KRONECKER
%
%  [vPhi, vE, vQ, vGam, vD, vK, vi] = simi(theta,din,thetau,dinu, T, N, i, ki)

i = i(:);
ix = size(i,1);

n = sum(ki);

[Phi, Gam, E, H, D, C, Q] = thd2ee(theta, din);
m = size(H,1);
r = max(size(Gam,2),size(D,2));

vPhi = zeros(N,n*m*ix);
vE   = zeros(N,n*m*ix);
vQ   = zeros(N,m*m*ix);
vK   = zeros(N,1);
vi   = zeros(N,ix);
vGam = zeros(N,n*r*ix);
vD   = zeros(N,m*r*ix);

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
             vK(j) = k;
             fmax = f;
          end

      if r
         [Phi, H, Tr, E, G] = echelon(ki,Phi, H, E, Gam, D);
         vGam(j,(k-1)*n*r+1:k*n*r) = G(:)';
         vD(j,(k-1)*m*r+1:k*m*r) = D(:)';
      else
         [Phi, H, Tr, E] = echelon(ki,Phi, H, E);
      end
      F = Phi(:,1:m);
      vPhi(j,(k-1)*n*m+1:k*n*m) = F(:)';
      vE(j,(k-1)*n*m+1:k*n*m)   = E(:)';
      vQ(j,(k-1)*m*m+1:k*m*m)   = Q(:)';
    end
end
sidopt('verbose','si');
