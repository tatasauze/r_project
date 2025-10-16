function  [n] = nida(z, nvec, s)


[N m] = size(z);
i = max(round(log(N)),max(nvec)+1);  
nv = size(nvec,2); 

% Estimating canonical correlations:  
S1 = singval(z, i, 0, s); % 2 = canonical, 1 = not canonical

%Ct = exp(-1.664)*N^(-.89)*i^(1.638);  % nid2
Ct = exp(.33)*N^(-.81424);
Pos = zeros(i,1);
   for k = 1:i
       dn = 2*(k-1)*m;     % Hannan-Deistler 
       Pos(k,1) = S1(k)^2 + dn*Ct;
   end
   
[~,a2]=min(Pos); n = a2-1;