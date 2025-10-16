function [T1, T2, vPhi, vE, vQ, porc] = tabvar(theta,din,lab,thetau,dinu,labu, T, N, i, ki)

disp('Modelo del output:');
prtmod(theta,din,lab)
if ~isempty(thetau)
   disp('Modelo del input:');
end

T = T(:);

fprintf(1,'Tamaños serie       : ');
fprintf(1,'%d ',T);
fprintf(1,'\nNúmero realizaciones: %d\n',N);
fprintf(1,'Valores de i        : %d al %d\n', i(1), i(max(size(i))));
fprintf(1,'Indices de Kronecker: ');
fprintf(1,'%d ',ki);
fprintf(1,'\n');

T1 = []; T2 = []; S1 = []; S2 = []; E1 = []; E2 = [];

porc = zeros(size(T));

for k=1:size(T,1)
    [vPhi,vsPhi,vH,vsH,vE,vsE,vQ,vGam,vsGam,vD,vsD, vK] = simvar(theta,din,thetau,dinu, T(k), N, i, ki);
    [mX1, sX1, eX1, porc]  = formod(theta,din,vPhi,vE,vQ,vGam,vD,vK,1,0);
    [mX2, sX2, eX2]  = formod(theta,din,vsPhi,vsE,zeros(size(vQ)),vsGam,vsD,vK,1,0);
    T1 = [T1;sX1*sqrt((N-1)/N)];
    T2 = [T2;mX2];
end