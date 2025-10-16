diary tablas.out
diary on

[tm1, dm1] = arma2thd([],[],[.5],[],[1],1);
[tm2, dm2] = arma2thd([],[],[-.5],[],[1],1);
[tm3, dm3] = arma2thd([],[],[-.9],[],[1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H, D, C, Q2] = thd2ee(tm2,dm2);
[Phi3, Gam, E3, H, D, C, Q3] = thd2ee(tm3,dm3);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);
E3 = E2-Phi3(:,1);

load SMA1100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
[mX3, sX3, eX3]  = formato(vPhim3,vEm3,vQm3,vKm3,Phi3,E3,Q3);
save TMA1100 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
disp('SMA1100')
mX1, sX1, eX1, mX2, sX2, eX2, mX3, sX3, eX3

load SMA1200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
[mX3, sX3, eX3]  = formato(vPhim3,vEm3,vQm3,vKm3,Phi3,E3,Q3);
save TMA1200 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
disp('SMA1200')
mX1, sX1, eX1, mX2, sX2, eX2, mX3, sX3, eX3

load SMA1000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
[mX3, sX3, eX3]  = formato(vPhim3,vEm3,vQm3,vKm3,Phi3,E3,Q3);
save TMA1000 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
disp('SMA1000')
mX1, sX1, eX1, mX2, sX2, eX2, mX3, sX3, eX3

[tm1, dm1] = arma2thd([],[],[-1.42 .73],[],[1],1);
[tm2, dm2] = arma2thd([],[],[-1.8 .9],[],[1],1);

[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H, D, C, Q2] = thd2ee(tm2,dm2);
Phi1 = Phi1(:,1);
Phi2 = Phi2(:,1);
E1 = E1-Phi1;
E2 = E2-Phi2;


load SMA2100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TMA2100 mX1 sX1 eX1 mX2 sX2 eX2
disp('SMA2100')
mX1, sX1, eX1, mX2, sX2, eX2

load SMA2200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TMA2200 mX1 sX1 eX1 mX2 sX2 eX2
disp('SMA2200')
mX1, sX1, eX1, mX2, sX2, eX2

load SMA2000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TMA2000 mX1 sX1 eX1 mX2 sX2 eX2
disp('SMA2000')
mX1, sX1, eX1, mX2, sX2, eX2

[tm1, dm1] = arma2thd([-.8],[],[-.5],[],[1],1);
[tm2, dm2] = arma2thd([-.5],[],[-.8],[],[1],1);
[tm3, dm3] = arma2thd([.5],[],[-.5],[],[1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H, D, C, Q2] = thd2ee(tm2,dm2);
[Phi3, Gam, E3, H, D, C, Q3] = thd2ee(tm3,dm3);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);
E3 = E2-Phi3(:,1);

load SARM1100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
[mX3, sX3, eX3]  = formato(vPhim3,vEm3,vQm3,vKm3,Phi3,E3,Q3);
save TARM1100 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
disp('SARM1100')
mX1, sX1, eX1, mX2, sX2, eX2, mX3, sX3, eX3

load SARM1200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
[mX3, sX3, eX3]  = formato(vPhim3,vEm3,vQm3,vKm3,Phi3,E3,Q3);
save TARM1200 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
disp('SARM1200')
mX1, sX1, eX1, mX2, sX2, eX2, mX3, sX3, eX3

load SARM1000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
[mX3, sX3, eX3]  = formato(vPhim3,vEm3,vQm3,vKm3,Phi3,E3,Q3);
save TARM1000 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
disp('SARM1000')
mX1, sX1, eX1, mX2, sX2, eX2, mX3, sX3, eX3

[tm1, dm1] = arma2thd([-.6],[],[-1 .64],[],[1],1);
[tm2, dm2] = arma2thd([-.6],[],[0 .64],[],[1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H, D, C, Q2] = thd2ee(tm2,dm2);
Phi1 = Phi1(:,1);
Phi2 = Phi2(:,1);
E1 = E1-Phi1;
E2 = E2-Phi2;

load SARM2100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TARM2100 mX1 sX1 eX1 mX2 sX2 eX2
disp('SARM2100')
mX1, sX1, eX1, mX2, sX2, eX2

load SARM2200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TARM2200 mX1 sX1 eX1 mX2 sX2 eX2
disp('SARM2200')
mX1, sX1, eX1, mX2, sX2, eX2

load SARM2000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TARM2000 mX1 sX1 eX1 mX2 sX2 eX2
disp('SARM2000')
mX1, sX1, eX1, mX2, sX2, eX2

[tm1, dm1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
[tm2, dm2] = arma2thd([0 .64],[],[-.6],[],[1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H, D, C, Q2] = thd2ee(tm2,dm2);
Phi1 = Phi1(:,1);
Phi2 = Phi2(:,1);
E1 = E1-Phi1;
E2 = E2-Phi2;

load SARM3100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TARM3100 mX1 sX1 eX1 mX2 sX2 eX2
disp('SARM3100')
mX1, sX1, eX1, mX2, sX2, eX2

load SARM3200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TARM3200 mX1 sX1 eX1 mX2 sX2 eX2
disp('SARM3200')
mX1, sX1, eX1, mX2, sX2, eX2

load SARM3000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
[mX2, sX2, eX2]  = formato(vPhim2,vEm2,vQm2,vKm2,Phi2,E2,Q2);
save TARM3000 mX1 sX1 eX1 mX2 sX2 eX2
disp('SARM3000')
mX1, sX1, eX1, mX2, sX2, eX2

[tm1, dm1] = arma2thd([-.7 0  0;0 0 0;0 -.4 0],[],[0 1.1 0; 0 -.6 0; 0 0 .5],[],[1 -.7 .4;-.7 1 0;.4 0 1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
E1 = E1-Phi1;

load SMVA100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
save TMVA100 mX1 sX1 eX1
disp('SMVA100')
mX1, sX1, eX1

load SMVA200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
save TMVA200 mX1 sX1 eX1
disp('SMVA200')
mX1, sX1, eX1

load SMVA000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
save TMVA000 mX1 sX1 eX1
disp('SMVA000')
mX1, sX1, eX1


[tm1, dm1] = arma2thd([-.4 -.3 .6;0 -.8 -.4;-.3 0 0],[],[-.7 0 0; -.1 -.2 0; .4 -.5 .1],[],[1 .5 .4;.5 1 .7;.4 .7 1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
E1 = E1-Phi1;

load SMVA2100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
save TMVA2100 mX1 sX1 eX1
disp('SMVA2100')
mX1, sX1, eX1

load SMVA2200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
save TMVA2200 mX1 sX1 eX1
disp('SMVA2200')
mX1, sX1, eX1

load SMVA2000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);
save TMVA2000 mX1 sX1 eX1
disp('SMVA2000')
mX1, sX1, eX1

diary off