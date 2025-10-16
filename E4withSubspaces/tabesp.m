function [T1, T2, porc, TPk1, TPk2] = tabesp(theta,din,lab,thetau,dinu,labu, T, N, i, ki, option, optval)

if nargin < 12, option = []; optval = 1; end
sopt = size(optval,1);

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
fprintf(1,'\nOpción de estim.  : %s\n', option);
fprintf(1,'Posibles valores    :');
for k=1:sopt
    fprintf(1,'%s, ',optval(k,:));
end
fprintf(1,'\n');

T1 = []; T2 = []; S1 = []; S2 = []; E1 = []; E2 = [];
TPk1 = []; TPk2 = []; SPk1 = []; SPk2 = []; EPk1 = []; EPk2 = [];
porc = zeros(size(T));

for k=1:size(T,1)
    if nargout > 3
       m = din(1,2);
       [vPhi, vE, vQ, vGam, vD, vK, vi,vPk] = simesp(theta,din,thetau,dinu, T(k), N, i, ki, option, optval);
       pk = size(vPk,2)/2;
       [mPk1, sPk1, ePk1]  = formod(theta,din,vPk(:,1:pk),-vPk(:,pk+1:2*pk),zeros(N,m*m),[],[],[],1,0);
       [mPk2, sPk2, ePk2]  = formod(theta,din,vPk(:,1:pk),-vPk(:,pk+1:2*pk),zeros(N,m*m),[],[],vK,1,0);
       TPk1 = [TPk1;mPk1]; SPk1 = [SPk1;sPk1]; EPk1 = [EPk1;ePk1];
       TPk2 = [TPk2;mPk2]; SPk2 = [SPk2;sPk2]; EPk2 = [EPk2;ePk2];
    else
       [vPhi, vE, vQ, vGam, vD, vK, vi] = simesp(theta,din,thetau,dinu, T(k), N, i, ki, option, optval);
    end
    [mX1, sX1, eX1]  = formod(theta,din,vPhi,vE,vQ,vGam,vD,[],sopt,0);
    [mX2, sX2, eX2, porc(k)]  = formod(theta,din,vPhi,vE,vQ,vGam,vD,vK,sopt,0);
    T1 = [T1;mX1]; S1 = [S1;sX1]; E1 = [E1;eX1];
    T2 = [T2;mX2]; S2 = [S2;sX2]; E2 = [E2;eX2];
end

T1 = [T1;S1;E1]; T2 = [T2;S2;E2];
TPk1 = [TPk1;SPk1;EPk1]; TPk2 = [TPk2;SPk2;EPk2];