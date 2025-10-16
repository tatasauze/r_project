load TERV1100
mT1 = mX1;
mT2 = mX2;
sT1 = sX1;
sT2 = sX2;
eT1 = eX1;
eT2 = eX2;
k11 = sk1;
k12 = sk2;

load TERV1200
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
sT1 = [sT1;sX1];
sT2 = [sT2;sX2];
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
k11 = [k11;sk1];
k12 = [k12;sk2];

load TERV1000
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
sT1 = [sT1;sX1]*sqrt(199/200);
sT2 = [sT2;sX2]*sqrt(199/200);
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
k11 = [k11;sk1];
k12 = [k12;sk2];

T1 = [mT1;sT1;eT1;mT2;sT2;eT2];

load TERV2100
mT1 = mX1;
mT2 = mX2;
sT1 = sX1;
sT2 = sX2;
eT1 = eX1;
eT2 = eX2;
k21 = sk1;
k22 = sk2;

load TERV2200
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
sT1 = [sT1;sX1];
sT2 = [sT2;sX2];
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
k21 = [k21;sk1];
k22 = [k22;sk2];

load TERV2000
mT1 = [mT1;mX1];
mT2 = [mT2;mX2];
sT1 = [sT1;sX1]*sqrt(199/200);
sT2 = [sT2;sX2]*sqrt(199/200);
eT1 = [eT1;eX1];
eT2 = [eT2;eX2];
k21 = [k21;sk1];
k22 = [k22;sk2];


T2 = [mT1;sT1;eT1;mT2;sT2;eT2];


