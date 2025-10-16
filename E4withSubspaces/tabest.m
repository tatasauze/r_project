
[tm1, dm1] = arma2thd([],[],[-.5],[-.5],[1],12);
[tm2, dm2] = arma2thd([],[],[-.7],[-.9],[1],12);
[tm3, dm3] = arma2thd([],[],[-.9],[-.7],[1],12);

load SSMA1100
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[0],[-.5],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.9],[0],[-.7],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.7],[0],[-.9],[1]);
save TSMA1100 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

load SSMA1200
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[0],[-.5],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.9],[0],[-.7],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.7],[0],[-.9],[1]);
save TSMA1200 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

load SSMA1000
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[0],[-.5],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.9],[0],[-.7],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.7],[0],[-.9],[1]);
save TSMA1000 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

[tm1, dm1] = arma2thd([-.5],[],[],[-.5],[1],12);
[tm2, dm2] = arma2thd([-.7],[],[],[-.9],[1],12);
[tm3, dm3] = arma2thd([-.9],[],[],[-.7],[1],12);

load SSAR1100
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[.5],[0],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.9],[.7],[0],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.7],[.9],[0],[1]);
save TSAR1100 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

load SSAR1200
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[.5],[0],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.9],[.7],[0],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.7],[.9],[0],[1]);
save TSAR1200 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

load SSAR1000
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[.5],[0],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.9],[.7],[0],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.7],[.9],[0],[1]);
save TSAR1000 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

[tm1, dm1] = arma2thd([-1 .64],[],[],[-.5],[1],12);
[tm2, dm2] = arma2thd([-1 .64],[],[],[-.7],[1],12);
[tm3, dm3] = arma2thd([-1 .64],[],[],[-.9],[1],12);

%[tm1, dm1] = arma2thd([-.3 .4],[],[],[-.5],[1],12);
%[tm2, dm2] = arma2thd([-.3 .4],[],[],[-.7],[1],12);
%[tm3, dm3] = arma2thd([-.3 .4],[],[],[-.9],[1],12);

load SSAR2100
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[1 -.64],[0 0],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.7],[1 -.64],[0 0],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.9],[1 -.64],[0 0],[1]);
save TSAR2100 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

load SSAR2200
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[1 -.64],[0 0],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.7],[1 -.64],[0 0],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.9],[1 -.64],[0 0],[1]);
save TSAR2200 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3

load SSAR2000
[mX1, sX1, eX1]  = format2(vPhis1,vEs1,vPhir1,vEr1,vQr1,[0],[-.5],[1 -.64],[0 0],[1]);
[mX2, sX2, eX2]  = format2(vPhis2,vEs2,vPhir2,vEr2,vQr2,[0],[-.7],[1 -.64],[0 0],[1]);
[mX3, sX3, eX3]  = format2(vPhis3,vEs3,vPhir3,vEr3,vQr3,[0],[-.9],[1 -.64],[0 0],[1]);
save TSAR2000 mX1 sX1 eX1 mX2 sX2 eX2 mX3 sX3 eX3
