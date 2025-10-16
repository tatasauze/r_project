
N = 200
sete4opt('econ','cero','vcon','cero');

[te de] = arma2thd([],[],[],[],[.5],1);
[tm1, dm1] = arma2thd([-1 .64],[],[],[],[1],1);
[tm1,dm1] = stackthd(tm1,dm1,te,de);
[tm1,dm1] = comp2thd(tm1,dm1);
[tm2, dm2] = arma2thd([0 .64],[],[],[],[1],1);
[tm2,dm2] = stackthd(tm2,dm2,te,de);
[tm2,dm2] = comp2thd(tm2,dm2);

sidinit
sidopt('pon','var');
sidopt('obs','re');
sidopt('met','exa');


[vPhi1, vE1, vQ11, vQ21, vR1] = simerv(tm1,dm1,100, N, [3:6], 2);
[vPhi2, vE2, vQ12, vQ22, vR2] = simerv(tm2,dm2,100, N, [3:6], 2);
save ERV1100 vPhi1 vE1 vPhi2 vE2 vQ11 vQ21 vQ12 vQ22 vR1 vR2

[vPhi1, vE1, vQ11, vQ21, vR1] = simerv(tm1,dm1,200, N, [3:6], 2);
[vPhi2, vE2, vQ12, vQ22, vR2] = simerv(tm2,dm2,200, N, [3:6], 2);
save ERV1200 vPhi1 vE1 vPhi2 vE2 vQ11 vQ21 vQ12 vQ22 vR1 vR2

[vPhi1, vE1, vQ11, vQ21, vR1] = simerv(tm1,dm1,1000, N, [3:6], 2);
[vPhi2, vE2, vQ12, vQ22, vR2] = simerv(tm2,dm2,1000, N, [3:6], 2);
save ERV1000 vPhi1 vE1 vPhi2 vE2 vQ11 vQ21 vQ12 vQ22 vR1 vR2

[te de] = arma2thd([],[],[],[],[.5],1);
[tm1, dm1] = arma2thd([-.5],[],[],[],[1],1);
[tm1,dm1] = comp2thd(tm1,dm1,[],[],te,de);
[tm2, dm2] = arma2thd([-.9],[],[],[],[1],1);
[tm2,dm2] = comp2thd(tm2,dm2,[],[],te,de);

[vPhi1, vE1, vQ11, vQ21, vR1] = simerv(tm1,dm1,100, N, [2:6], 1);
[vPhi2, vE2, vQ12, vQ22, vR2] = simerv(tm2,dm2,100, N, [2:6], 1);
save ERV2100 vPhi1 vE1 vPhi2 vE2 vQ11 vQ21 vQ12 vQ22 vR1 vR2

[vPhi1, vE1, vQ11, vQ21, vR1] = simerv(tm1,dm1,200, N, [2:6], 1);
[vPhi2, vE2, vQ12, vQ22, vR2] = simerv(tm2,dm2,200, N, [2:6], 1);
save ERV2200 vPhi1 vE1 vPhi2 vE2 vQ11 vQ21 vQ12 vQ22 vR1 vR2

[vPhi1, vE1, vQ11, vQ21, vR1] = simerv(tm1,dm1,1000, N, [2:6], 1);
[vPhi2, vE2, vQ12, vQ22, vR2] = simerv(tm2,dm2,1000, N, [2:6], 1);
save ERV2000 vPhi1 vE1 vPhi2 vE2 vQ11 vQ21 vQ12 vQ22 vR1 vR2
