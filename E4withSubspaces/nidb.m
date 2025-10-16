function  [n] = nidb(z, nvec)

[N m] = size(z);
i = max(round(log(N)),max(nvec)+1); 
nv = size(nvec,2); 

S1 = singval(z, i, 1, 1); % 2 = canonical, 1 = not canonical

Ct = exp(-1.995)*N^(-.87)*i^(1.63);  % nid2
Pos = zeros(i,1);
   for k = 1:i
       dn = 2*(k-1)*m;     % Hannan-Deistler 
       Pos(k,1) = S1(k)^2 + dn*Ct;
   end
   
[~,a2] = min(Pos); n = a2-1;