
N = 200
sete4opt('econ','cero','vcon','cero');

sidinit

[tu, du] = arma2thd([-.7],[],[],[],[1],1);
[tm1, dm1] = arma2thd([],[],[-.7],[],[1],1, [.5], 1);
[tm2, dm2] = arma2thd([],[],[-.7],[],[1],1, [.4 .5], 1);

load EX1100
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 100, N, [2;3;4;5;10], 1);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 100, N, [3;4;5;6;10], 1);
save EX1100 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

load EX1200
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 200, N, [2;3;4;5;10], 1);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 200, N, [3;4;5;6;10], 1);
save EX1200 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

load EX1000
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 1000, N, [2;3;4;5;10], 1);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 1000, N, [3;4;5;6;10], 1);
save EX1000 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

[tm1, dm1] = arma2thd([-1 .64],[],[-.6],[],[1],1, [.5],1);
[tm2, dm2] = arma2thd([-1 .64],[],[-.6],[],[1],1, [0 0 .5],1);

load EX2100
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 100, N, [3;4;5;6;10], 2);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 100, N, [4;5;6;7;10], 2);
save EX2100 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

load EX2200
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 200, N, [3;4;5;6;10], 2);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 200, N, [4;5;6;7;10], 2);
save EX2200 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

load EX2000
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 1000, N, [3;4;5;6;10], 2);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 1000, N, [4;5;6;7;10], 2);
save EX2000 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

[tm1, dm1] = arma2thd([0 .64],[],[-.6],[],[1],1,[.5],1); 
[tm2, dm2] = arma2thd([0 .64],[],[-.6],[],[1],1,[.4 0 .5],1);

load EX3100
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 100, N, [3;4;5;6;10], 2);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 100, N, [4;5;6;7;10], 2);
save EX3100 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

load EX3200
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 200, N, [3;4;5;6;10], 2);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 200, N, [4;5;6;7;10], 2);
save EX3200 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

load EX3000
%[vPhim1, vEm1, vQm1, vGm1, vDm1, vKm1] = simtoex(tm1,dm1, tu, du, 1000, N, [3;4;5;6;10], 2);
[vPhim2, vEm2, vQm2, vGm2, vDm2, vKm2] = simtoex(tm2,dm2, tu, du, 1000, N, [4;5;6;7;10], 2);
save EX3000 vPhim1 vEm1 vQm1 vGm1 vDm1 vKm1 vPhim2 vEm2 vQm2 vGm2 vDm2 vKm2

