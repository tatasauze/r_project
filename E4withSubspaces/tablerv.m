[te de] = arma2thd([],[],[],[],[.5],1);
[tm1, dm1] = arma2thd([-1 .64],[],[],[],[1],1);
[tm1,dm1] = comp2thd(tm1,dm1,[],[],te,de);
[tm2, dm2] = arma2thd([0 .64],[],[],[],[1],1);
[tm2,dm2] = comp2thd(tm2,dm2,[],[],te,de);

[Phi1, Gam, E1, H1, D, C, Q21, S1, R1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H2, D, C, Q22, S2, R2] = thd2ee(tm2,dm2);
R1 = .5; R2 = .5;
P1 = lyapunov(Phi1,E1*Q21*E1');
G1 = Phi1*P1*H1' + E1*Q21;
L01 = H1*P1*H1' + 1.5;
[P1, E1, Q11] = Ricnit(Phi1, H1, G1, L01);
P2 = lyapunov(Phi2,E2*Q22*E2');
G2 = Phi2*P2*H2' + E2*Q22;
L02 = H2*P2*H2' + 1.5;
[P2, E2, Q12] = Ricnit(Phi2, H2, G2, L02);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);


load ERV1100
k1 = find(vR1 >= 0); k2 = find(vR2 >= 0);
sk1 = size(k1,1)/200; sk2 = size(k2,1)/200;
[mX1, sX1, eX1]  = formater(vPhi1(k1,:),vE1(k1,:),vQ11(k1,:),vQ21(k1,:),vR1(k1,:),Phi1,E1,Q11,Q21,R1);
[mX2, sX2, eX2]  = formater(vPhi2(k2,:),vE2(k2,:),vQ12(k2,:),vQ22(k2,:),vR2(k2,:),Phi2,E2,Q12,Q22,R2);
save TERV1100 mX1 sX1 eX1 mX2 sX2 eX2  sk1 sk2

load ERV1200
k1 = find(vR1 >= 0); k2 = find(vR2 >= 0);
sk1 = size(k1,1)/200; sk2 = size(k2,1)/200;
[mX1, sX1, eX1]  = formater(vPhi1(k1,:),vE1(k1,:),vQ11(k1,:),vQ21(k1,:),vR1(k1,:),Phi1,E1,Q11,Q21,R1);
[mX2, sX2, eX2]  = formater(vPhi2(k2,:),vE2(k2,:),vQ12(k2,:),vQ22(k2,:),vR2(k2,:),Phi2,E2,Q12,Q22,R2);
save TERV1200 mX1 sX1 eX1 mX2 sX2 eX2 sk1 sk2

load ERV1000
k1 = find(vR1 >= 0); k2 = find(vR2 >= 0);
sk1 = size(k1,1)/200; sk2 = size(k2,1)/200;
[mX1, sX1, eX1]  = formater(vPhi1(k1,:),vE1(k1,:),vQ11(k1,:),vQ21(k1,:),vR1(k1,:),Phi1,E1,Q11,Q21,R1);
[mX2, sX2, eX2]  = formater(vPhi2(k2,:),vE2(k2,:),vQ12(k2,:),vQ22(k2,:),vR2(k2,:),Phi2,E2,Q12,Q22,R2);
save TERV1000 mX1 sX1 eX1 mX2 sX2 eX2 sk1 sk2

[te de] = arma2thd([],[],[],[],[.5],1);
[tm1, dm1] = arma2thd([-.5],[],[],[],[1],1);
[tm1,dm1] = comp2thd(tm1,dm1,[],[],te,de);
[tm2, dm2] = arma2thd([-.9],[],[],[],[1],1);
[tm2,dm2] = comp2thd(tm2,dm2,[],[],te,de);

[Phi1, Gam, E1, H1, D, C, Q21, S1, R1] = thd2ee(tm1,dm1);
[Phi2, Gam, E2, H2, D, C, Q22, S2, R2] = thd2ee(tm2,dm2);
R1 = .5; R2 = .5;
P1 = lyapunov(Phi1,E1*Q21*E1');
G1 = Phi1*P1*H1' + E1*Q21;
L01 = H1*P1*H1' + 1.5;
[P1, E1, Q11] = Ricnit(Phi1, H1, G1, L01);
P2 = lyapunov(Phi2,E2*Q22*E2');
G2 = Phi2*P2*H2' + E2*Q22;
L02 = H2*P2*H2' + 1.5;
[P2, E2, Q12] = Ricnit(Phi2, H2, G2, L02);
E1 = E1-Phi1(:,1);
E2 = E2-Phi2(:,1);


load ERV2100
k1 = find(vR1 >= 0); k2 = find(vR2 >= 0);
sk1 = size(k1,1)/200; sk2 = size(k2,1)/200;
[mX1, sX1, eX1]  = formater(vPhi1(k1,:),vE1(k1,:),vQ11(k1,:),vQ21(k1,:),vR1(k1,:),Phi1,E1,Q11,Q21,R1);
[mX2, sX2, eX2]  = formater(vPhi2(k2,:),vE2(k2,:),vQ12(k2,:),vQ22(k2,:),vR2(k2,:),Phi2,E2,Q12,Q22,R2);
save TERV2100 mX1 sX1 eX1 mX2 sX2 eX2 sk1 sk2

load ERV2200
k1 = find(vR1 >= 0); k2 = find(vR2 >= 0);
sk1 = size(k1,1)/200; sk2 = size(k2,1)/200;
[mX1, sX1, eX1]  = formater(vPhi1(k1,:),vE1(k1,:),vQ11(k1,:),vQ21(k1,:),vR1(k1,:),Phi1,E1,Q11,Q21,R1);
[mX2, sX2, eX2]  = formater(vPhi2(k2,:),vE2(k2,:),vQ12(k2,:),vQ22(k2,:),vR2(k2,:),Phi2,E2,Q12,Q22,R2);
save TERV2200 mX1 sX1 eX1 mX2 sX2 eX2 sk1 sk2

load ERV2000
k1 = find(vR1 >= 0); k2 = find(vR2 >= 0);
sk1 = size(k1,1)/200; sk2 = size(k2,1)/200;
[mX1, sX1, eX1]  = formater(vPhi1(k1,:),vE1(k1,:),vQ11(k1,:),vQ21(k1,:),vR1(k1,:),Phi1,E1,Q11,Q21,R1);
[mX2, sX2, eX2]  = formater(vPhi2(k2,:),vE2(k2,:),vQ12(k2,:),vQ22(k2,:),vR2(k2,:),Phi2,E2,Q12,Q22,R2);
save TERV2000 mX1 sX1 eX1 mX2 sX2 eX2 sk1 sk2
