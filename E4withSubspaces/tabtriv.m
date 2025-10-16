[tm1, dm1] = arma2thd([-.4 -.3 .6;0 -.8 -.4;-.3 0 0],[],[-.7 0 0; -.1 -.2 0; .4 -.5 .1],[],[1 .5 .4;.5 1 .7;.4 .7 1],1);
[Phi1, Gam, E1, H, D, C, Q1] = thd2ee(tm1,dm1);
E1 = E1-Phi1;

load SMVA2100
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);

T100 = [];
E100 = [];
S100 = [];

k = [1:4:36 37:8:180 41:8:180];

for j=0:3
    for i=1:5
        T = reshape(mX1(i,k+j),3,size(mX1,2)/12);
        S = reshape(sX1(i,k+j),3,size(sX1,2)/12);
        E = reshape(eX1(i,k+j),3,size(eX1,2)/12);
        T100 = [T100;T];
        S100 = [S100;S];
        E100 = [E100;E];
    end
end

load SMVA2200
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);

T200 = [];
E200 = [];
S200 = [];

for j=0:3
    for i=1:5
        T = reshape(mX1(i,k+j),3,size(mX1,2)/12);
        S = reshape(sX1(i,k+j),3,size(sX1,2)/12);
        E = reshape(eX1(i,k+j),3,size(eX1,2)/12);
        T200 = [T200;T];
        S200 = [S200;S];
        E200 = [E200;E];
    end
end


load SMVA2000
[mX1, sX1, eX1]  = formato(vPhim1,vEm1,vQm1,vKm1,Phi1,E1,Q1);

T000 = [];
E000 = [];
S000 = [];

for j=0:3
    for i=1:5
        T = reshape(mX1(i,k+j),3,size(mX1,2)/12);
        S = reshape(sX1(i,k+j),3,size(sX1,2)/12);
        E = reshape(eX1(i,k+j),3,size(eX1,2)/12);
        T000 = [T000;T];
        S000 = [S000;S];
        E000 = [E000;E];
    end
end
