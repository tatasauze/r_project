% Analisis de número de incendios
% TF con inputs: Precio de la Cebada (pts/Hl)
%                Denuncias
%                Lluvia media
%                Has. de monte público
% Modelo compuesto
clear
clc
eeinit;

load incendi.da1;
% Precio cebada
x1=incendi(:,9);
x1(1:34,1)=x1(1:34,1)-mean(x1(1:34,1));
% Denuncias totales
x2=incendi(:,6)./100;
x2(4:37,1)=x2(4:37,1)-mean(x2(4:37,1));
% Precipitación media
x3=incendi(:,5)./100;
x3=x3-mean(x3);
% Has. de monte público
x4=incendi(:,10);
x4=x4-mean([x4(3:18);x4(20:35);x4(37:37)]);
% Número de incendios
y=incendi(:,2);
y(1:36,1)=y(1:36,1)-mean(y(1:36,1));

% Modelo del output
[t1,d1,l1] = tf2thd([-.64 .24 -.38],[],[],[],[2.41],1, ...
[1.0 NaN;1.81 NaN;.50 NaN;NaN 1.0],[]);
% Modelo de los inputs
phi1=[-.81 NaN NaN NaN; NaN NaN NaN NaN; NaN NaN NaN NaN;NaN NaN NaN NaN];
phi2=[ .74 NaN NaN NaN; NaN NaN NaN NaN; NaN NaN NaN NaN;NaN NaN NaN NaN];
phi3=[-.53 NaN NaN NaN; NaN NaN NaN NaN; NaN NaN NaN NaN;NaN NaN NaN NaN];
sigma=[1.06 .06 -.11 .0; .06 1.0 -.14 .0; -.11 -.14 1.91 .0; 0 0 0 5.0];
[t2,d2,l2] = arma2thd([phi1 phi2 phi3],[],[],[],[sigma],1);
% Modelo compuesto
[theta,din] = comp2thd(t1,d1,t2,d2);
lab=[l1;l2];

sete4opt('vcond','lyap','econd','ml','var','fac','maxit',500);
theta=e4preest(theta,din,[y x1 x2 x3 x4]);
prtmod(theta,din,lab);

tic
[thopt,it,fval,g,h]=eemin('fvmiss',theta,'',din,[y x1 x2 x3 x4]);
prtest(thopt, din, lab, it, fval, g, h, [], [],(toc/60));
