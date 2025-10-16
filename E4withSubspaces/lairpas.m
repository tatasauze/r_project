% AIRPAS.M  Serie de pasajeros de avión Box y Jenkins 1976
% 
% Tamaño de la serie: 144
% diferencias: 1 regular y 1 estacional (s=12)
%
% especificación aceptada: media móvil doble
%
% valores atípicos: 29 (I), 54 (E), 135 (I) y 62 (E, hallada en 2ª vuelta)

!del airpas.res
sidinit
load airpas.ser
y = airpas;

diary airpas.res
diary on

y = transdif(y,0,1,1,12);

[tchi2s, desvs, angles, rvchis, dgs] = sidang(y,[],3,12)

[Phis,sPhis,H,sH,Es,sEs,Q,yr] = sidechel(y, [], [2], 1, 12);
yr = yr';
Phis, sPhis
Es, sEs

[tchi2, desv, angle, rvchi, dg] = sidang(yr,[],6)
[Phi1,sPhi1,H,sH,E1,sE1,Q1,innov1] = sidechei(yr, [], [2:4], 1);
Phi1, sPhi1, E1, sE1, Q1
[tchia1, desva1, anglea1, rvchia1, dga1] = sidang(innov1',[],6)
[tchias1, desvas1, angleas1, rvchias1, dgas1] = sidang(innov1',[],3,12)
disp('***************************************************')
int = zeros(144,5);
int(29,1) = 1;
int(39:144,2)  = ones(106,1);
int(54:144,3) = ones(91,1);
int(62:144,4)  = ones(83,1);
int(135,5) = 1;

u = transdif(int,1,1,1,12);

[Phis,sPhis,H,sH,Es,sEs,Q,Ds,sDs,yr] = sidint(y, u, [2], 1, 12);
yr = yr';
Phis, sPhis
Es, sEs
Ds, sDs

[tchi2, desv, angle, rvchi, dg] = sidang(yr,[],6)
[Phi2,sPhi2,H,sH,E2,sE2,Q2,innov2] = sidechei(yr, [], [2:4], 1);

Phi2, sPhi2, E2, sE2, Q2

[tchia2, desva2, anglea2, rvchia2, dga2] = sidang(innov2',[],6)
[tchias2, desvas2, angleas2, rvchias2, dgas2] = sidang(innov2',[],3,12)

uidents(innov1');
uidents(innov2');
plotsers(innov1');
plotsers(innov2');

diary off