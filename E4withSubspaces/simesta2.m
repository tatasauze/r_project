function [vPhis, vEs, vPhir, vEr, vQr, vK] = simestat(theta,din,T, N, is, ir, kis, kir, s)
%
%  Simula modelos estacionales con matrices block-Hankel extendidas y no extendidas
%  Atención, fuera de la función es preciso poner las opciones de estimación por 
%  defecto
%
% [vPhi, vE, vQ, vK] = simtot(theta,din,T, N, i, n, s, ext)

if nargin < 9,   s = 1; end
is = is(:); ixs = size(is,1);
ir = ir(:); ixr = size(ir,1);

ns = sum(kis);
nr = sum(kir);

[Phi,ign,E,H] = thd2ee(theta,din);
m = size(H,1);

vPhis = zeros(N,ns*m*2);
size(vPhis)
vEs   = zeros(N,ns*m*2);
vK    = zeros(N,3);

vPhir = zeros(N,nr*m*2);
vEr   = zeros(N,nr*m*2);
vQr   = zeros(N,m*m*2);

ext = ['n';'s'];
sidopt('verbose','no');

for j=1:N
%
j
    y = simmod(theta,din,T);
    fmax = inf;
    ptrs = 0;
    ptrr = 0;
    ptr2 = 0;
    sidopt('ext',ext(1,:));
    t = sidang(y, [], is(ixs),s);
    k = find(diag(t) < .95);
    if isempty(k), k = is(ixs)*m; else k = k(1)-1; end
    vK(j,1) = k;
    if k == ns
       for k1=1:2
           sidopt('ext',ext(k1,:));
           for k=1:size(is,1)
               [Phis,Hs,Es,Qs,innov] = sident(y, [], is(k), ns, s);
               t = sidang(innov', [], ir(ixr));
               kix = find(diag(t) < .95);
               if isempty(kix), kix = ir(ixr)*m; else kix = kix(1)-1; end
               if kix == nr
                  vK(j,k1+1) = kix;
                  for k2=1:size(ir,1)
                    [Phir,Hr,Er,Qr,innov2] = sident(innov', [], ir(k2), nr);
                    [thet, di] = ee2thd(Phir, [], Er, Hr, [], [], Qr);
                    f = pem(thet,di,innov');
                    if f < fmax
                       fmax = f;
                       Phisx = Phis;
                       Hsx   = Hs;
                       Esx   = Es;
                       Phirx = Phir;
                       Hrx   = Hr;
                       Erx   = Er;
                       Qrx   = Qr; 
                    end
                  end
               end
           end
           if vK(j,k1+1) == nr
              [Phi, H, ign, E] = echelon(kis,Phisx, Hsx, Esx);
              F = Phi(:,1:m);
              vPhis(j, ptrs+1:ptrs+ns*m) = F(:)';
              vEs(j, ptrs+1:ptrs+ns*m)   = E(:)';
              [Phi, H, ign, E] = echelon(kir,Phirx, Hrx, Erx);
              F = Phi(:,1:m);
              vPhir(j, ptrr+1:ptrr+nr*m) = F(:)';
              vEr(j, ptrr+1:ptrr+nr*m)   = E(:)';
              vQr(j, ptr2+1:ptr2+m*m)   = Qrx(:)';
           end
           ptrs = ptrs +ns*m;
           ptrr = ptrr +nr*m;
           ptr2 = ptr2 + m*m;
       end
    end
end
sidopt('verbose','si');
