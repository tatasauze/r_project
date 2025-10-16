function [T1, T2] = tabint(theta,din,lab,u, T, N, i, ki, option, optval)

if nargin < 12, option = []; optval = 1; end
sopt = size(optval,1);

disp('Modelo del output:');
prtmod(theta,din,lab)

T = T(:);

fprintf(1,'Tamaños serie       : ');
fprintf(1,'%d ',T);
fprintf(1,'\nNúmero realizaciones: %d\n',N);
fprintf(1,'Valores de i        : %d al %d\n', i(1), i(max(size(i))));
fprintf(1,'Indices de Kronecker: ');
fprintf(1,'%d ',ki);
fprintf(1,'\nOpción de estim.  : %s\n', option);
fprintf(1,'Posibles valores    :');
for k=1:sopt
    fprintf(1,'%s, ',optval(k,:));
end
fprintf(1,'\n');

T1 = []; T2 = []; S1 = []; S2 = []; E1 = []; E2 = [];
porc = zeros(size(T));

for k=1:size(T,1)

    [vPhi, vsPhi, vE, vsE, vQ, vD, vsD, vD2] = simint(theta,din,u(1:T(k),:),T(k), N, i, ki);
    [mX1, sX1, eX1]  = formint(theta,din,vPhi,vE,vQ,vD,vD2,[],sopt,0);
    mX2 = formint(theta,din,vsPhi,vsE,vQ,vsD,vsD,[],sopt,0);

    T1 = [T1;mX1]; S1 = [S1;sX1]; E1 = [E1;eX1];
    T2 = [T2;mX2];
end

T1 = [T1;S1;E1]; 
