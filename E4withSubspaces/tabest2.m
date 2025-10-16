load tsma1100
mT1 = mX1;
mT2 = mX2;
mT3 = mX3;
sT1 = sX1;
sT2 = sX2;
sT3 = sX3;
eT1 = eX1;
eT2 = eX2;
eT3 = eX3;

load tsma1200
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
mT3 = [mT3;mX3];
sT1 = [sT1;sX1];
sT2 = [sT2;sX2];
sT3 = [sT3;sX3];
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
eT3 = [eT3;eX3];

load tsma1000
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
mT3 = [mT3;mX3];
sT1 = [sT1;sX1]*sqrt(199/200);
sT2 = [sT2;sX2]*sqrt(199/200);
sT3 = [sT3;sX3]*sqrt(199/200);
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
eT3 = [eT3;eX3];


T1 = [mT1;sT1;eT1;mT2;sT2;eT2;mT3;sT3;eT3];

load tsar1100
mT1 = mX1;
mT2 = mX2;
mT3 = mX3;
sT1 = sX1;
sT2 = sX2;
sT3 = sX3;
eT1 = eX1;
eT2 = eX2;
eT3 = eX3;

load tsar1200
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
mT3 = [mT3;mX3];
sT1 = [sT1;sX1];
sT2 = [sT2;sX2];
sT3 = [sT3;sX3];
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
eT3 = [eT3;eX3];

load tsar1000
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
mT3 = [mT3;mX3];
sT1 = [sT1;sX1]*sqrt(199/200);
sT2 = [sT2;sX2]*sqrt(199/200);
sT3 = [sT3;sX3]*sqrt(199/200);
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
eT3 = [eT3;eX3];

T2 = [mT1;sT1;eT1;mT2;sT2;eT2;mT3;sT3;eT3];

load tsar2100
mT1 = mX1;
mT2 = mX2;
mT3 = mX3;
sT1 = sX1;
sT2 = sX2;
sT3 = sX3;
eT1 = eX1;
eT2 = eX2;
eT3 = eX3;

load tsar2200
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
mT3 = [mT3;mX3];
sT1 = [sT1;sX1];
sT2 = [sT2;sX2];
sT3 = [sT3;sX3];
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
eT3 = [eT3;eX3];

load tsar2000
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
mT3 = [mT3;mX3];
sT1 = [sT1;sX1]*sqrt(199/200);
sT2 = [sT2;sX2]*sqrt(199/200);
sT3 = [sT3;sX3]*sqrt(199/200);
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
eT3 = [eT3;eX3];

T3 = [mT1;sT1;eT1;mT2;sT2;eT2;mT3;sT3;eT3];