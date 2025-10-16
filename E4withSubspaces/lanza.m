
N = 200
sete4opt('econ','cero','vcon','cero');

%[tm1, dm1] = arma2thd([],[],[.5],[],[1],1);
%[tm2, dm2] = arma2thd([],[],[-.5],[],[1],1);
%[tm3, dm3] = arma2thd([],[],[-.9],[],[1],1);

sidinit

% load SMA1100
%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [2;3;4;5;10;15], 1);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,100, N, [2;3;4;5;10;15], 1);
%[vPhim3, vEm3, vQm3, vKm3] = simtot(tm3,dm3,100, N, [2;3;4;5;10;15], 1);

%save SMA1100 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2 vPhim3 vEm3 vQm3 vKm3

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [2;3;4;5;10;15], 1);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,200, N, [2;3;4;5;10;15], 1);
%[vPhim3, vEm3, vQm3, vKm3] = simtot(tm3,dm3,200, N, [2;3;4;5;10;15], 1);

%save SMA1200 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2 vPhim3 vEm3 vQm3 vKm3

%load SMA1000
%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [2;3;4;5;10;15], 1);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,1000, N, [2;3;4;5;10;15], 1);
%[vPhim3, vEm3, vQm3, vKm3] = simtot(tm3,dm3,1000, N, [2;3;4;5;10;15], 1);

%save SMA1000 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2 vPhim3 vEm3 vQm3 vKm3

%[tm1, dm1] = arma2thd([],[],[-1.42 .73],[],[1],1);
%[tm2, dm2] = arma2thd([],[],[-1.8 .9],[],[1],1);

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,100, N, [3;4;5;6;10;15], 2);

%save SMA2100 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,200, N, [3;4;5;6;10;15], 2);

%save SMA2200 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,1000, N, [3;4;5;6;10;15], 2);

%save SMA2000 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[tm1, dm1] = arma2thd([-.8],[],[-.5],[],[1],1);
%[tm2, dm2] = arma2thd([-.5],[],[-.8],[],[1],1);
%[tm3, dm3] = arma2thd([.5],[],[-.5],[],[1],1);

%sidinit

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [2;3;4;5;10;15], 1);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,100, N, [2;3;4;5;10;15], 1);
%[vPhim3, vEm3, vQm3, vKm3] = simtot(tm3,dm3,100, N, [2;3;4;5;10;15], 1);

%save SARM1100 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2 vPhim3 vEm3 vQm3 vKm3

%load SARM1200
%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [2;3;4;5;10;15], 1);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,200, N, [2;3;4;5;10;15], 1);
%[vPhim3, vEm3, vQm3, vKm3] = simtot(tm3,dm3,200, N, [2;3;4;5;10;15], 1);

%save SARM1200 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2 vPhim3 vEm3 vQm3 vKm3

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [2;3;4;5;10;15], 1);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,1000, N, [2;3;4;5;10;15], 1);
%[vPhim3, vEm3, vQm3, vKm3] = simtot(tm3,dm3,1000, N, [2;3;4;5;10;15], 1);

%save SARM1000 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2 vPhim3 vEm3 vQm3 vKm3

%[tm1, dm1] = arma2thd([-.6],[],[-1 .64],[],[1],1);
%[tm2, dm2] = arma2thd([-.6],[],[0 .64],[],[1],1);

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,100, N, [3;4;5;6;10;15], 2);

%save SARM2100 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,200, N, [3;4;5;6;10;15], 2);

%save SARM2200 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,1000, N, [3;4;5;6;10;15], 2);

%save SARM2000 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[tm1, dm1] = arma2thd([-1 .64],[],[-.6],[],[1],1);
%[tm2, dm2] = arma2thd([0 .64],[],[-.6],[],[1],1);

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,100, N, [3;4;5;6;10;15], 2);

%save SARM3100 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,200, N, [3;4;5;6;10;15], 2);

%save SARM3200 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [3;4;5;6;10;15], 2);
%[vPhim2, vEm2, vQm2, vKm2] = simtot(tm2,dm2,1000, N, [3;4;5;6;10;15], 2);

%save SARM3000 vPhim1 vEm1 vQm1 vKm1 vPhim2 vEm2 vQm2 vKm2

%[tm1, dm1] = arma2thd([[-.7 0  0;0 0 0;0 -.4 0],[],[0 1.1 0; 0 -.6 0; 0 0 .5],[],[1 -.7 .4;-.7 1 0;.4 0 1],1);

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [3;4;5;6;10], [1;1;1]);

%save SMVA100 vPhim1 vEm1 vQm1 vKm1

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [3;4;5;6;10], [1;1;1]);

%save SMVA200 vPhim1 vEm1 vQm1 vKm1

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [3;4;5;6;10], [1;1;1]);

%save SMVA000 vPhim1 vEm1 vQm1 vKm1

disp('hola')

[tm1, dm1] = arma2thd([-.4 -.3 .6;0 -.8 -.4;-.3 0 0],[],[-.7 0 0; -.1 -.2 0; .4 -.5 .1],[],[1 .5 .4;.5 1 .7;.4 .7 1],1);

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,100, N, [3;4;5;6;10], [1;1;1]);

%save SMVA2100 vPhim1 vEm1 vQm1 vKm1

%[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,200, N, [3;4;5;6;10], [1;1;1]);

%save SMVA2200 vPhim1 vEm1 vQm1 vKm1

[vPhim1, vEm1, vQm1, vKm1] = simtot(tm1,dm1,1000, N, [3;4;5;6;10], [1;1;1]);

save SMVA2000 vPhim1 vEm1 vQm1 vKm1