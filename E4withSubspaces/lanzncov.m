% Harvey (1989)

Nr = 200;
Phi = [1 1;0 1];
E = eye(2);
H = [1 0];
C = 1;
Gam = [];
D = [];
Q = diag([2;1]);
R = [.5];
S = [];
[Betah1, VBetah1] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 100, Nr, 3);
[Betah2, VBetah2] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 200, Nr, 3);
[Betah5, VBetah5] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 500, Nr, 3);
[Betah10, VBetah10] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 1000, Nr, 3);
k = sqrt((Nr-1)/Nr);
tablah = [mean(Betah1); mean(Betah2); mean(Betah5); mean(Betah10);std(Betah1)*k; std(Betah2)*k; std(Betah5)*k; std(Betah10)*k ;mean(VBetah1); mean(VBetah2); mean(VBetah5); mean(VBetah10)];
save simncovh tablah Betah1 Betah2 Betah5 Betah10 VBetah1 VBetah2 VBetah5 VBetah10 

% Young(1988)

Phi = [1 1;0 1];
E = [0;1];
H = [1 0];
C = 1;
Gam = [];
D = [];
Q = 2;
R = .5;
S = [];
[Betay1, VBetay1] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 100, Nr, 2);
[Betay2, VBetay2] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 200, Nr, 2);
[Betay5, VBetay5] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 500, Nr, 2);
[Betay10, VBetay10] = simncov(Phi, Gam, E, H, D, C, Q, S, R, 34, 1000, Nr, 2);
tablay = [mean(Betay1); mean(Betay2); mean(Betay5); mean(Betay10);std(Betay1)*k; std(Betay2)*k; std(Betay5)*k; std(Betay10)*k ;mean(VBetay1); mean(VBetay2); mean(VBetay5); mean(VBetay10)];
save simncovy tablay Betay1 Betay2 Betay5 Betay10 VBetay1 VBetay2 VBetay5 VBetay10 
