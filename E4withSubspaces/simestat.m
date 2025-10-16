function [vPhis, vEs, vPhir, vEr, vQr, vK] = simestat(theta,din,T, N, is, ir, kis, kir, s)
%
%  Simula modelos estacionales con matrices block-Hankel extendidas y no extendidas
%  Atención, fuera de la función es preciso poner las opciones de estimación por 
%  defecto
%
% [vPhi, vE, vQ, vK] = simtot(theta,din,T, N, i, n, s, ext)

if nargin < 8, ext = 0; end
if nargin < 7,   s = 1; end
is = is(:);
ir = ir(:);

ns = sum(kis);
nr = sum(kir);

[Phi,ign,E,H] = thd2ee(theta,din);
m = size(H,1);

vPhis = zeros(N,ns*m);
vEs   = zeros(N,ns*m);
vQs   = zeros(N,m*m);
vK    = zeros(N,is*2+1);

vPhir = zeros(N,nr*m);
vEr   = zeros(N,nr*m);
vQr   = zeros(N,m*m);

ext = ['n';'s'];
sidopt('verbose','no');

for j=1:N
%
    y = simmod(theta,din,T);
    ptrs = 0;
    ptrr = 0;
    ptr2 = 0;
    t = sidang(y, [], is(size(is,1)),12);
    kix = find(diag(t) < .95);

    if isempty(kix), kix = is(size(is,1))*m; else kix = kix(1)-1; end
    vK(j,1) = kix;

    for k1=1:2
        fmax = inf;
        sidopt('ext',ext(k1,:));
        for k=1:size(is,1)
            [Phis,Hs,Es,Qs,innov] = sident(y, [], is(k), ns, s);

            t = sidang(innov', [], ir(size(ir,1)));
            kix = find(diag(t) < .95);
            if isempty(kix), kix = ir(size(ir,1))*m; else kix = kix(1)-1; end
            vK(j,size(is,1)*(k1-1)+k+1) = kix;

            for k2=1:size(ir,1)
                [Phir,Hr,Er,Qr,innov2] = sident(innov', [], ir(k2), nr);
                [thet, di] = ee2thd(Phir, [], Er, Hr, [], [], Qr);
                f = pem(thet,di,innov');
                if f < fmax
                   fmax = f;
%               if trace(innov2*innov2'/size(innov,2)) < fmax
%                  fmax = trace(innov2*innov2'/size(innov,2));
                  Phisx = Phis;
                  Hsx   = Hs;
                  Esx   = Es;

%        Qsx = Qs;
                  Phirx = Phir;
                  Hrx   = Hr;
                  Erx   = Er;
                  Qrx   = Qr; 
               end
            end
        end
        [Phi, H, ign, E] = echelon(kis,Phisx, Hsx, Esx);
        F = Phi(:,1:m);
        vPhis(j, ptrs+1:ptrs+ns*m) = F(:)';
        vEs(j, ptrs+1:ptrs+ns*m)   = E(:)';

%        vQs(j, ptr2+1:ptr2+m*m)   = Qsx(:)';

        [Phi, H, ign, E] = echelon(kir,Phirx, Hrx, Erx);
        F = Phi(:,1:m);
        vPhir(j, ptrr+1:ptrr+nr*m) = F(:)';
        vEr(j, ptrr+1:ptrr+nr*m)   = E(:)';
        vQr(j, ptr2+1:ptr2+m*m)   = Qrx(:)';
        ptrs = ptrs +ns*m;
        ptrr = ptrr +nr*m;
        ptr2 = ptr2 + m*m;
    end
end
sidopt('verbose','si');
