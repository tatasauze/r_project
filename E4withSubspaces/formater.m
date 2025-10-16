function [mX1, sX1, eX1]  = format2(vPhi,vE,vQ1,vQ2,vR,Phi,E,Q1,Q2,R);
%

   Phi = Phi(:,1)';
   E   = E(:)';
   Phie = sqrt(mean((vPhi - kron(ones(size(vPhi,1),1),Phi)).^2));
   Ee   = sqrt(mean((vE - kron(ones(size(vE,1),1),E)).^2));
   Q1e  = sqrt(mean((vQ1 - kron(ones(size(vQ1,1),1),Q1)).^2));
   Q2e  = sqrt(mean((vQ2 - kron(ones(size(vQ2,1),1),Q2)).^2));
   Re   = sqrt(mean((vR  - kron(ones(size(vR,1),1),R)).^2));

   Phim = mean(vPhi);
   Em   = mean(vE);
   Q1m   = mean(vQ1);
   Q2m   = mean(vQ2);
   Rm   = mean(vR);

   Phis = std(vPhi);
   Es   = std(vE);
   Q1s  = std(vQ1);
   Q2s  = std(vQ2);
   Rs   = std(vR);


   mX1 = [Phim Em Q1m Q2m Rm];
   sX1 = [Phis Es Q1s Q2s Rs];
   eX1 = [Phie Ee Q1e Q2e Re];
