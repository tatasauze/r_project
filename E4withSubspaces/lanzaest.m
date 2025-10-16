
N = 200
sete4opt('econ','cero','vcon','cero');

[tm1, dm1] = arma2thd([],[],[-.5],[-.5],[1],12);
[tm2, dm2] = arma2thd([],[],[-.7],[-.9],[1],12);
[tm3, dm3] = arma2thd([],[],[-.9],[-.7],[1],12);

sidinit
sidopt('pon','var');
sidopt('obs','re');
sidopt('met','exa');

%[vPhis1, vEs1, vPhir1, vEr1, vQr1,vK1] = simestat(tm1,dm1,100, N, [2;3], [2:6], 1, 1, 12);
%[vPhis2, vEs2, vPhir2, vEr2, vQr2,vK2] = simestat(tm2,dm2,100, N, [2;3], [2:6], 1, 1, 12);
%[vPhis3, vEs3, vPhir3, vEr3, vQr3,vK3] = simestat(tm3,dm3,100, N, [2;3], [2:6], 1, 1, 12);
%save SSMA1100 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

%[vPhis1, vEs1, vPhir1, vEr1, vQr1,vK1] = simestat(tm1,dm1,200, N, [2;3], [2:6], 1, 1, 12);
%[vPhis2, vEs2, vPhir2, vEr2, vQr2,vK2] = simestat(tm2,dm2,200, N, [2;3], [2:6], 1, 1, 12);
%[vPhis3, vEs3, vPhir3, vEr3, vQr3,vK3] = simestat(tm3,dm3,200, N, [2;3], [2:6], 1, 1, 12);
%save SSMA1200 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

%[vPhis1, vEs1, vPhir1, vEr1, vQr1,vK1] = simestat(tm1,dm1,1000, N, [2;3], [2:6], 1, 1, 12);
%[vPhis2, vEs2, vPhir2, vEr2, vQr2,vK2] = simestat(tm2,dm2,1000, N, [2;3], [2:6], 1, 1, 12);
%[vPhis3, vEs3, vPhir3, vEr3, vQr3,vK3] = simestat(tm3,dm3,1000, N, [2;3], [2:6], 1, 1, 12);
%save SSMA1000 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

[tm1, dm1] = arma2thd([-.5],[],[],[-.5],[1],12);
[tm2, dm2] = arma2thd([-.7],[],[],[-.9],[1],12);
[tm3, dm3] = arma2thd([-.9],[],[],[-.7],[1],12);

load SSAR1100
%[vPhis1, vEs1, vPhir1, vEr1, vQr1,vK1] = simestat(tm1,dm1,100, N, [2;3], [2:6], 1, 1, 12);
%[vPhis2, vEs2, vPhir2, vEr2, vQr2,vK2] = simestat(tm2,dm2,100, N, [2;3], [2:6], 1, 1, 12);
[vPhis3, vEs3, vPhir3, vEr3, vQr3,vK3] = simestat(tm3,dm3,100, N, [2;3], [2:6], 1, 1, 12);
save SSAR1100 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

load SSAR1200
%[vPhis1, vEs1, vPhir1, vEr1, vQr1, vK1] = simestat(tm1,dm1,200, N, [2;3], [2:6], 1, 1, 12);
%[vPhis2, vEs2, vPhir2, vEr2, vQr2, vK2] = simestat(tm2,dm2,200, N, [2;3], [2:6], 1, 1, 12);
[vPhis3, vEs3, vPhir3, vEr3, vQr3, vK3] = simestat(tm3,dm3,200, N, [2;3], [2:6], 1, 1, 12);
save SSAR1200 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

load SSAR1000
%[vPhis1, vEs1, vPhir1, vEr1, vQr1, vK1] = simestat(tm1,dm1,1000, N, [2;3], [2:6], 1, 1, 12);
%[vPhis2, vEs2, vPhir2, vEr2, vQr2, vK2] = simestat(tm2,dm2,1000, N, [2;3], [2:6], 1, 1, 12);
[vPhis3, vEs3, vPhir3, vEr3, vQr3, vK3] = simestat(tm3,dm3,1000, N, [2;3], [2:6], 1, 1, 12);
save SSAR1000 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3


return

[tm1, dm1] = arma2thd([-1 .64],[],[],[-.5],[1],12);
[tm2, dm2] = arma2thd([-1 .64],[],[],[-.7],[1],12);
[tm3, dm3] = arma2thd([-1 .64],[],[],[-.9],[1],12);

[vPhis1, vEs1, vPhir1, vEr1, vQr1, vK1] = simestat(tm1,dm1,100, N, [2;3], [3:6], 1, 2, 12);
[vPhis2, vEs2, vPhir2, vEr2, vQr2, vK2] = simestat(tm2,dm2,100, N, [2;3], [3:6], 1, 2, 12);
[vPhis3, vEs3, vPhir3, vEr3, vQr3, vK3] = simestat(tm3,dm3,100, N, [2;3], [3:6], 1, 2, 12);
save SSAR2100 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

[vPhis1, vEs1, vPhir1, vEr1, vQr1, vK1] = simestat(tm1,dm1,200, N, [2;3], [3:6], 1, 2, 12);
[vPhis2, vEs2, vPhir2, vEr2, vQr2, vK2] = simestat(tm2,dm2,200, N, [2;3], [3:6], 1, 2, 12);
[vPhis3, vEs3, vPhir3, vEr3, vQr3, vK3] = simestat(tm3,dm3,200, N, [2;3], [3:6], 1, 2, 12);
save SSAR2200 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

[vPhis1, vEs1, vPhir1, vEr1, vQr1, vK1] = simestat(tm1,dm1,1000, N, [2;3], [3:6], 1, 2, 12);
[vPhis2, vEs2, vPhir2, vEr2, vQr2, vK2] = simestat(tm2,dm2,1000, N, [2;3], [3:6], 1, 2, 12);
[vPhis3, vEs3, vPhir3, vEr3, vQr3, vK3] = simestat(tm3,dm3,1000, N, [2;3], [3:6], 1, 2, 12);
save SSAR2000 vPhis1 vEs1 vPhir1 vEr1 vQr1 vPhis2 vEs2 vPhir2 vEr2 vQr2 vPhis3 vEs3 vPhir3 vEr3 vQr3 vK1 vK2 vK3

