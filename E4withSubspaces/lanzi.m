% Simulación total utilizando simesp
diary LANZi
diary on

sidinit
N = 200
sete4opt('econ','cero','vcon','cero');
sidopt('met','exac','obs','reestim','ext','si', 'pon','v');

[tu, du,lu] = arma2thd([-.7],[],[],[],[1],1);
[tm1, dm1,l1] = arma2thd([],[],[.5],[],[1],1);
[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1);
[tm3, dm3,l3] = arma2thd([],[],[-.9],[],[1],1);
T1 = tabi(tm1,dm1,l1,[],[],[], [100;200;1000], N, [2:6 8 10], 1)
T2 = tabi(tm2,dm2,l2,[],[],[], [100;200;1000], N, [2:6 8 10], 1)
T3 = tabi(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6 8 10], 1)
save iMA1 T1 T2 T2 T3

diary off
diary on

[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6],[],[1],1);
T1 = tabi(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6 8 10], 2)
T2 = tabi(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6 8 10], 2)
save iAR2 T1 T2

diary off
diary on

[tm1, dm1,l1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);
T1 = tabi(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [2:6 8 10], 1)
T2 = tabi(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [2:6 8 10], 1)
save iMAX1 T1 T2
diary off
diary on

[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
[tm2, dm2,l2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);
T1 = tabi(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [3:6 8 10], 2)
T2 = tabi(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [3:6 8 10], 2)
save iARX2 T1 T2

diary off
