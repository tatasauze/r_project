% Simulación total utilizando simesp
diary LANZAVAR
diary on

sidinit
N = 1000
sete4opt('econ','cero','vcon','cero');
sidopt('met','exac','obs','reestim','ext','si');

[tu, du,lu] = arma2thd([-.7],[],[],[],[1],1);
[tm1, dm1,l1] = arma2thd([],[],[.5],[],[1],1);
[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1);
[tm3, dm3,l3] = arma2thd([],[],[-.9],[],[1],1);
[M1, V1, vPhi1, vE1, vQ1, porc1] = tabvar(tm1,dm1,l1,[],[],[], [100;200;1000], N, [2:6], 1);
[M2, V2, vPhi2, vE2, vQ2, porc2] = tabvar(tm2,dm2,l2,[],[],[], [100;200;1000], N, [2:6], 1);
[M3, V3, vPhi3, vE3, vQ3, porc3] = tabvar(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6], 1);
save VAR_MA1 M1 M2 M3 V1 V2 V3 vPhi1 vE1 vQ1 vPhi2 vE2 vQ2 vPhi3 vE3 vQ3 porc1 porc2 porc3

diary off
diary on

[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6],[],[1],1);
[M1, V1, vPhi1, vE1, vQ1, porc1] = tabvar(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], 2);
[M2, V2, vPhi2, vE2, vQ2, porc2] = tabvar(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6], 2);
save VAR_AR2 M1 M2 V1 V2 vPhi1 vE1 vQ1 vPhi2 vE2 vQ2 porc1 porc2

diary off
diary on

[tm1, dm1,l1] = arma2thd([-.4 -.3 .6;0 -.8 -.4;-.3 0 0],[],[-.7 0 0; -.1 -.2 0; .4 -.5 .1],[],[1 .5 .4;.5 1 .7;.4 .7 1],1);
[M1, V1, vPhi1, vE1, vQ1, porc] = tabvar(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], [1;1;1]);
save VAR_VMA M1 V1 vPhi1 vE1 vQ1 porc
diary off
diary on

