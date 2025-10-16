function [vPhi,vsPhi,vH,vsH,vE,vsE,vQ,vGam,vsGam,vD,vsD,vK] = simvar(theta,din,thetau,dinu, T, N, i, ki)
%
%  Simulación matrices de varianzas
%

i = i(:);
ix = size(i,1);

n = sum(ki);

[Phi, Gam, E, H, D, C, Q] = thd2ee(theta, din);
m = size(H,1);
r = max(size(Gam,2),size(D,2));

mb = sum(ki > 0);

vPhi = zeros(N,n*mb);
vE   = zeros(N,n*m);
vH   = zeros(N,m*mb);
vsPhi = zeros(N,n*mb);
vsE   = zeros(N,n*m);
vsH   = zeros(N,m*mb);
vQ   = zeros(N,m*m);
if r
   vGam = zeros(N,n*r);
   vD   = zeros(N,m*r);
   vsGam = zeros(N,n*r);
   vsD   = zeros(N,m*r);
end
vK = zeros(N,1);

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
    
    if r
    else
       [Phi,sPhi,H,sH,E,sE,Q] = sidechei(y, [], i, ki);
    end

    F = Phi(:,1:mb);
    vPhi(j,:) = F(:)';
    vE(j,:)   = E(:)';
    F = sPhi(:,1:mb);
    vsPhi(j,:) = F(:)';
    vsE(j,:)   = sE(:)';
    vQ(j,:)   = Q(:)';
end
sidopt('verbose','si');
