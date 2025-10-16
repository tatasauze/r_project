% Simulación total utilizando simesp
%diary LANZATOT
%diary on

sidinit
N = 200
sete4opt('econ','cero','vcon','cero');
sidopt('met','exac','obs','reestim','ext','si');

[tu, du,lu] = arma2thd([-.7],[],[],[],[1],1);
%[tm1, dm1,l1] = arma2thd([],[],[.5],[],[1],1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1);
%[tm3, dm3,l3] = arma2thd([],[],[-.9],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [2:6], 1, 'pon', ['var';'no '])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [2:6], 1, 'pon', ['var';'no '])
%[T13, T23, porc3] = tabesp(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6], 1, 'pon', ['var';'no '])
%save MA1PON T11 T21 porc1 T12 T22 porc2 T13 T23 porc3

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
%[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], 2, 'pon', ['var';'no '])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6], 2, 'pon', ['var';'no '])
%save AR2PON T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [2:6], 1, 'pon', ['var';'no '])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [2:6], 1, 'pon', ['var';'no '])
%save MAX1PON T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
%[tm2, dm2,l2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [3:6], 2, 'pon', ['var';'no '])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [3:6], 2, 'pon', ['var';'no '])
%save ARX2PON T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%sidopt('met','exac','obs','reestim','pon','v');

%[tm1, dm1,l1] = arma2thd([],[],[.5],[],[1],1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1);
%[tm3, dm3,l3] = arma2thd([],[],[-.9],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [2:6], 1, 'ext', ['si';'no'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [2:6], 1, 'ext', ['si';'no'])
%[T13, T23, porc3] = tabesp(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6], 1, 'ext', ['si';'no'])
%save MA1EXT T11 T21 porc1 T12 T22 porc2 T13 T23 porc3

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
%[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], 2, 'ext', ['si';'no'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6], 2, 'ext', ['si';'no'])
%save AR2EXT T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [2:6], 1, 'ext', ['si';'no'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [2:6], 1, 'ext', ['si';'no'])
%save MAX1EXT T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
%[tm2, dm2,l2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [3:6], 2, 'ext', ['si';'no'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [3:6], 2, 'ext', ['si';'no'])
%save ARX2EXT T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%sidopt('met','exac','ext','si','pon','v');

%[tm1, dm1,l1] = arma2thd([],[],[.5],[],[1],1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1);
%[tm3, dm3,l3] = arma2thd([],[],[-.9],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [2:6], 1, 'obs', ['reest';'moore'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [2:6], 1, 'obs', ['reest';'moore'])
%[T13, T23, porc3] = tabesp(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6], 1, 'obs', ['reest';'moore'])
%save MA1OBS T11 T21 porc1 T12 T22 porc2 T13 T23 porc3

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
%[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], 2, 'obs', ['reest';'moore'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6], 2, 'obs', ['reest';'moore'])
%save AR2OBS T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [2:6], 1, 'obs', ['reest';'moore'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [2:6], 1, 'obs', ['reest';'moore'])
%save MAX1OBS T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
%[tm2, dm2,l2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [3:6], 2, 'obs', ['reest';'moore'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [3:6], 2, 'obs', ['reest';'moore'])
%save ARX2OBS T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%sidopt('obs','reest','ext','si','pon','v');

%[tm1, dm1,l1] = arma2thd([],[],[.5],[],[1],1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1);
%[tm3, dm3,l3] = arma2thd([],[],[-.9],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [2:6], 1, 'met', ['exact';'aprox'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [2:6], 1, 'met', ['exact';'aprox'])
%[T13, T23, porc3] = tabesp(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6], 1, 'met', ['exact';'aprox'])
%save MA1MET T11 T21 porc1 T12 T22 porc2 T13 T23 porc3

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
%[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6],[],[1],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], 2, 'met', ['exact';'aprox'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6], 2, 'met', ['exact';'aprox'])
%save AR2MET T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
%[tm2, dm2,l2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [2:6], 1, 'met', ['exact';'aprox'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [2:6], 1, 'met', ['exact';'aprox'])
%save MAX1MET T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

%[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
%[tm2, dm2,l2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);
%[T11, T21, porc1] = tabesp(tm1,dm1,l1,tu,du,lu, [100;200;1000], N, [3:6], 2, 'met', ['exact';'aprox'])
%[T12, T22, porc2] = tabesp(tm2,dm2,l2,tu,du,lu, [100;200;1000], N, [3:6], 2, 'met', ['exact';'aprox'])
%save ARX2MET T11 T21 porc1 T12 T22 porc2

%diary off
%diary on

load VARMAMET
[tm1, dm1,l1] = arma2thd([-.7 0  0;0 0 0;0 -.4 0],[],[0 1.1 0; 0 -.6 0; 0 0 .5],[],[1 -.7 .4;-.7 1 0;.4 0 1],1);
[tm2, dm2,l2] = arma2thd([-.4 -.3 .6;0 -.8 -.4;-.3 0 0],[],[-.7 0 0; -.1 -.2 0; .4 -.5 .1],[],[1 .5 .4;.5 1 .7;.4 .7 1],1);
ma1 = [-.8 -.7; NaN -.5];
[tm3, dm3, l3] = arma2thd([],[],ma1,[],[1.0 1.0],1);
ma1 = [-1.6 -.7; NaN -.5]; ma2 = [.63 NaN; NaN -.4];
[tm4, dm4, l4] = arma2thd([],[],[ma1 ma2],[],[1.0 1.0],1);
%[T11, T21, porc1,vPK11, vPK21] = tabesp(tm1,dm1,l1,[],[],[], [100;200;1000], N, [3:6], [1;1;1], 'met', ['exact';'aprox']);
%[T12, T22, porc2,vPK12, vPK22] = tabesp(tm2,dm2,l2,[],[],[], [100;200;1000], N, [3:6], [1;1;1], 'met', ['exact';'aprox']);
%[T13, T23, porc3,vPK13, vPK23] = tabesp(tm3,dm3,l3,[],[],[], [100;200;1000], N, [2:6], [1;1], 'met', ['exact';'aprox']);
[T14, T24, porc4,vPK14, vPK24] = tabesp(tm4,dm4,l4,[],[],[], [100;200;1000], N, [4:6], [2;2], 'met', ['exact';'aprox']);
save VARMAMET T11 T21 porc1 T12 T22 porc2 vPK11 vPK21  vPK12 vPK22 T13 T23 porc3 T14 T24 porc4 vPK13 vPK23  vPK14 vPK24 

%diary off
