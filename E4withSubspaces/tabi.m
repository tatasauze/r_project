function T  = tabi(theta,din,lab,thetau,dinu,labu, T, N, i, ki)

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
fprintf(1,'\n');

T1 = []; T2 = []; S1 = []; S2 = []; E1 = []; E2 = [];

porc = zeros(size(T));

n = sum(ki);
m = din(1,2);
r = din(1,3);
vPhix = zeros(N,n*m);
vEx   = zeros(N,n*m);
vQx   = zeros(N,m*m);
vGamx = zeros(N,n*r);
vDx   = zeros(N,m*r);

i = i(:);

for k=1:size(T,1)
    [vPhi, vE, vQ, vGam, vD, vK] = simi(theta,din,thetau,dinu, T(k), N, i, ki);
    [mX1, sX1, eX1]  = formod(theta,din,vPhi,vE,vQ,vGam,vD,[],size(i,1),0);
    for j=1:N
        vPhix(j,:) = vPhi(j,(vK(j)-1)*n*m+1:vK(j)*n*m);
        vEx(j,:) =   vE(j,(vK(j)-1)*n*m+1:vK(j)*n*m);
        vQx(j,:) =   vQ(j,(vK(j)-1)*m*m+1:vK(j)*m*m);
        if r
           vGamx(j,:) = vGam(j,(vK(j)-1)*n*r+1:vK(j)*n*r);
           vDx(j,:) =     vD(j,(vK(j)-1)*m*r+1:vK(j)*m*r);
        end
    end
    
    mX1 = reshape(mX1,size(mX1,2)/size(i,1),size(i,1))';
    sX1 = reshape(sX1,size(sX1,2)/size(i,1),size(i,1))';
    eX1 = reshape(eX1,size(eX1,2)/size(i,1),size(i,1))';

    [mX2, sX2, eX2]  = formod(theta,din,vPhix,vEx,vQx,vGamx,vDx,[],1,0);

    mX1 = [mX1;mX2];
    sX1 = [sX1;sX2];
    eX1 = [eX1;eX2];

    T1 = [T1;mX1]; S1 = [S1;sX1]; E1 = [E1;eX1];
end

T = [T1;S1;E1];