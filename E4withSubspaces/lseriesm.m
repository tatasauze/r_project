load seriesm.ser
sidinit
sidopt('ext','n');
!del seriesm.res
diary seriesm.res
diary on
y = transdif(seriesm,1,1);
u = y(:,1);
u = u-mean(u);
y = y(:,2);
y = y-mean(y);

[tchi2, desv, angle, rvchi, dg] = sidang(y,u,6)
[Phi,sPhi,H,sH,E,sE,Q,Gam,sGam,D,sD,innov] = sidechei(y, u, [4:8], 3);
[tchi2r, desvr, angler, rvchir, dgr] = sidang(innov' - mean(innov'),[],6)
Phi
sPhi
E
sE
Q
G=Gam-Phi(:,1)*D;
G
D

[Phi2,sPhi2,H2,sH2,E2,sE2,Q2,innov2] = sidechei(innov', [], [2:8], 1);
[tchi2r2, desvr2, angler2, rvchir2, dgr2] = sidang(innov2',[],6)
Phi2
sPhi2
E2
sE2
Q2
diary off