[tm1, dm1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
[tm2, dm2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);

[Phi1, Gam1, E1, H1, D1, C1, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam2, E2, H2, D2, C2, Q2] = thd2ee(tm2,dm2);
Phi1 = Phi1(:,1);
Phi2 = Phi2(:,1);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);
Gam1 = Gam1-Phi1(:,1)*D1;
Gam2 = Gam2-Phi2(:,1)*D2;

load EX1100
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX1100 mX1 sX1 eX1 mX2 sX2 eX2

load EX1200
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX1200 mX1 sX1 eX1 mX2 sX2 eX2

load EX1000
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX1000 mX1 sX1 eX1 mX2 sX2 eX2

[tm1, dm1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
[tm2, dm2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);

[Phi1, Gam1, E1, H1, D1, C1, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam2, E2, H2, D2, C2, Q2] = thd2ee(tm2,dm2);
Phi1 = Phi1(:,1);
Phi2 = Phi2(:,1);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);
Gam1 = Gam1-Phi1(:,1)*D1;
Gam2 = Gam2-Phi2(:,1)*D2;

load EX2100
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX2100 mX1 sX1 eX1 mX2 sX2 eX2

load EX2200
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX2200 mX1 sX1 eX1 mX2 sX2 eX2

load EX2000
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX2000 mX1 sX1 eX1 mX2 sX2 eX2

[tm1, dm1] = arma2thd([0 .64],[],[-.6],[],[1],1,[.5],1); 
[tm2, dm2] = arma2thd([0 .64],[],[-.6],[],[1],1,[.4 0 .5],1);


[Phi1, Gam1, E1, H1, D1, C1, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam2, E2, H2, D2, C2, Q2] = thd2ee(tm2,dm2);
Phi1 = Phi1(:,1);
Phi2 = Phi2(:,1);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);
Gam1 = Gam1-Phi1(:,1)*D1;
Gam2 = Gam2-Phi2(:,1)*D2;

load EX3100
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX3100 mX1 sX1 eX1 mX2 sX2 eX2

load EX3200
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX3200 mX1 sX1 eX1 mX2 sX2 eX2

load EX3000
[mX1, sX1, eX1]  = formatx(vPhim1,vEm1,vQm1,vDm1,vGm1,Phi1,E1,Q1,D1,Gam1);
[mX2, sX2, eX2]  = formatx(vPhim2,vEm2,vQm2,vDm2,vGm2,Phi2,E2,Q2,D2,Gam2);
save TEX3000 mX1 sX1 eX1 mX2 sX2 eX2

