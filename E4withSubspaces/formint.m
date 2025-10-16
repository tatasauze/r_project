function  [mX, sX, eX, porc]  = formint(theta,din,vPhi,vE,vQ,vD,vD2,vK,nopt,flag)
%
%  nopt = número de opciones
%  flag = 0 => tablas separadas
%  flag = 1 => tablas juntas

   [Phi, Gam, E, H, D, C, Q] = thd2ee(theta, din); 
   m = size(H,1);
   r = max(size(Gam,2),size(D,2));
   n = size(Phi,1);

   Phi = Phi(:,1:m);
   E = E-Phi;

   if isempty(vK)
      N = size(vPhi,1);
   else
      k = vK == n;
      if size(vK,2) > 1
         k = cumprod(k')';
      end
      k = find(k(:,size(k,2)) == 1) ;
      N = size(k,1);
      porc = size(k,1)/size(vPhi,1);
      vPhi = vPhi(k,:);
      vE   = vE(k,:);
      vQ   = vQ(k,:);
      if r
         vD = vD(k,:);
         vD2= vD(k,:);
      end
   end
   Phie = sqrt(mean((vPhi - kron(ones(N,nopt),Phi(:)')).^2));
   Ee   = sqrt(mean((vE   - kron(ones(N,nopt),E(:)')).^2));
   Qe   = sqrt(mean((vQ   - kron(ones(N,nopt),Q(:)')).^2));
   Phim = mean(vPhi);
   Phis = std(vPhi);
   Em   = mean(vE);
   Es   = std(vE);
   Qm   = mean(vQ);
   Qs   = std(vQ);
   if r
      Ge   = sqrt(mean((vD - kron(ones(N,nopt),D(:)')).^2));
      De   = sqrt(mean((vD2 - kron(ones(N,nopt),D(:)')).^2));
      Gm   = mean(vD);
      Gs   = std(vD);
      Dm   = mean(vD2);
      Ds   = std(vD2);
   end

   if m > 1
      Phie = blocrot(reshape(Phie,n,m*nopt)',m);
      Phim = blocrot(reshape(Phim,n,m*nopt)',m);
      Phis = blocrot(reshape(Phis,n,m*nopt)',m);
      Ee = blocrot(reshape(Ee,n,m*nopt)',m);
      Em = blocrot(reshape(Em,n,m*nopt)',m);
      Es = blocrot(reshape(Es,n,m*nopt)',m);
      Qe = blocrot(reshape(Qe,m,m*nopt)',m);
      Qm = blocrot(reshape(Qm,m,m*nopt)',m);
      Qs = blocrot(reshape(Qs,m,m*nopt)',m);
      if r
         Ge = blocrot(reshape(Ge,m,r*nopt)',r);
         Gm = blocrot(reshape(Gm,m,r*nopt)',r);
         Gs = blocrot(reshape(Gs,m,r*nopt)',r);
         De = blocrot(reshape(De,m,r*nopt)',r);
         Dm = blocrot(reshape(Dm,m,r*nopt)',r);
         Ds = blocrot(reshape(Ds,m,r*nopt)',r);
      end
   end
            
   if flag
      kn = kron(ones(1,nopt),[1:nopt:n*nopt]);
      kn = kn + kron([0:nopt-1],ones(1,n));
      km = kron(ones(1,nopt),[1:nopt:m*nopt]);
      km = km + kron([0:nopt-1],ones(1,m));
      Phim = Phim(:,kn);      
      Phie = Phie(:,kn);      
      Phis = Phis(:,kn);      
      Em = Em(:,kn);      
      Ee = Ee(:,kn);      
      Es = Es(:,kn);      
      Qm = Qm(:,km);      
      Qe = Qe(:,km);      
      Qs = Qs(:,km);      
      if r
         Dm = Dm(:,km);      
         De = De(:,km);      
         Ds = Ds(:,km);      
         Gm = Gm(:,km);      
         Ge = Ge(:,km);      
         Gs = Gs(:,km);      
         mX = [Phim Gm Dm Em Qm];
         sX = [Phis Gs Ds Es Qs];
         eX = [Phie Ge De Ee Qe];
      else
         mX = [Phim Em Qm];
         sX = [Phis Es Qs];
         eX = [Phie Ee Qe];
      end
   else
      mX = []; sX = []; eX = [];
      for i=1:nopt
          if r
             mX = [mX Phim(:,(i-1)*n+1:i*n) Gm(:,(i-1)*r+1:i*r) Dm(:,(i-1)*r+1:i*r) Em(:,(i-1)*n+1:i*n) Qm(:,(i-1)*m+1:i*m)];
             sX = [sX Phis(:,(i-1)*n+1:i*n) Gs(:,(i-1)*r+1:i*r) Ds(:,(i-1)*r+1:i*r) Es(:,(i-1)*n+1:i*n) Qs(:,(i-1)*m+1:i*m)];
             eX = [eX Phie(:,(i-1)*n+1:i*n) Ge(:,(i-1)*r+1:i*r) De(:,(i-1)*r+1:i*r) Ee(:,(i-1)*n+1:i*n) Qe(:,(i-1)*m+1:i*m)];
          else
             mX = [mX Phim(:,(i-1)*n+1:i*n) Em(:,(i-1)*m+1:i*m) Qm(:,(i-1)*m+1:i*m)];
             sX = [sX Phis(:,(i-1)*n+1:i*n) Es(:,(i-1)*m+1:i*m) Qs(:,(i-1)*m+1:i*m)];
             eX = [eX Phie(:,(i-1)*n+1:i*n) Ee(:,(i-1)*m+1:i*m) Qe(:,(i-1)*m+1:i*m)];
          end
      end
   end

