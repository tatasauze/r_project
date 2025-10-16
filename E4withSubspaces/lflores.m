% LFLORES.M  Precios de harina de Tiao y Tsay
% 
% Tamaño de la serie: 100x3
% diferencias: 1 regular
%

!del flores.res
diary flores.res
diary on

sidinit;
load flores.ser
y = flores;

y = transdif(y,0,1);

[tchi211, desv11, angle11, rvchi11, dg11] = sidang(y(:,1),[],4)
[tchi212, desv12, angle12, rvchi12, dg12] = sidang(y(:,2),[],4)
[tchi213, desv13, angle13, rvchi13, dg13] = sidang(y(:,3),[],4)

[tchi21, desv1, angle1, rvchi1, dg1] = sidang(y,[],3)

[Phi,sPhi,H,sH,E,sE,Q,yr] = sidechei(y, [], [2:5], [1;0;0]);
Phi, sPhi, H, sH, E, sE, Q
[tchi21r, desv1r, angle1r, rvchi1r, dg1r] = sidang(yr',[],3)

[Phi2,sPhi2,H2,sH2,E2,sE2,Q2,yr2] = sidechei(y, [], [2:5], [1;1;0]);
Phi2, sPhi2, H2, sH2, E2, sE2, Q2
[tchi22r, desv2r, angle2r, rvchi2r, dg2r] = sidang(yr2',[],3)

%uidents(yr');
%uidents(yr2');
%plotsers(yr');
%plotsers(yr2');

diary off