load TEX1100
T111 = [mX1;sX1;eX1];
T121 = [mX2;sX2;eX2];

s = sprintf('l2c1:l%1dc%1d',size(T111,1)+1,size(T111,2))
ddepoke(chan,s,T111)
j = size(T111,1)+2+1;
j2 = j+size(T121,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T121,2))
ddepoke(chan,s,T121)
j = j2+2;

load TEX1200
T112 = [mX1;sX1;eX1];
T122 = [mX2;sX2;eX2];

j2 = j+size(T112,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T112,2))
ddepoke(chan,s,T112)
j = j2+2;
j2 = j+size(T122,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T122,2))
ddepoke(chan,s,T122)
j = j2+2;

load TEX1000
T110 = [mX1;sX1;eX1];
T120 = [mX2;sX2;eX2];

j2 = j+size(T110,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T110,2))
ddepoke(chan,s,T110)
j = j2+2;
j2 = j+size(T120,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T120,2))
ddepoke(chan,s,T120)
j = j2+2;

load TEX2100
T211 = [mX1;sX1;eX1];
T221 = [mX2;sX2;eX2];

j2 = j+size(T211,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T211,2))
ddepoke(chan,s,T211)
j = j2+2;
j2 = j+size(T221,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T221,2))
ddepoke(chan,s,T221)
j = j2+2;

load TEX2200
T212 = [mX1;sX1;eX1];
T222 = [mX2;sX2;eX2];

j2 = j+size(T212,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T212,2))
ddepoke(chan,s,T212)
j = j2+2;
j2 = j+size(T222,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T222,2))
ddepoke(chan,s,T222)
j = j2+2;

load TEX2000
T210 = [mX1;sX1;eX1];
T220 = [mX2;sX2;eX2];

j2 = j+size(T210,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T210,2))
ddepoke(chan,s,T210)
j = j2+2;
j2 = j+size(T220,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T220,2))
ddepoke(chan,s,T220)
j = j2+2;

load TEX3100
T311 = [mX1;sX1;eX1];
T321 = [mX2;sX2;eX2];

j2 = j+size(T311,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T311,2))
ddepoke(chan,s,T311)
j = j2+2;
j2 = j+size(T321,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T321,2))
ddepoke(chan,s,T321)
j = j2+2;

load TEX3200
T312 = [mX1;sX1;eX1];
T322 = [mX2;sX2;eX2];

j2 = j+size(T312,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T312,2))
ddepoke(chan,s,T312)
j = j2+2;
j2 = j+size(T322,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T322,2))
ddepoke(chan,s,T322)
j = j2+2;

load TEX3000
T310 = [mX1;sX1;eX1];
T320 = [mX2;sX2;eX2];

j2 = j+size(T310,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T310,2))
ddepoke(chan,s,T310)
j = j2+2;
j2 = j+size(T320,1)-1;
s = sprintf('l%1dc1:l%1dc%1d',j,j2,size(T320,2))
ddepoke(chan,s,T320)
j = j2+2;


