function [mX1, sX1, eX1]  = format2(vPhis,vEs,vPhir,vEr,vQ,Phis,Es,Phir,Er,Q);
%

   Phies = sqrt(mean((vPhis - kron(ones(size(vPhis,1),2),Phis)).^2));
   Ees   = sqrt(mean((vEs - kron(ones(size(vEs,1),2),Es)).^2));
   Phier = sqrt(mean((vPhir - kron(ones(size(vPhir,1),2),Phir)).^2));
   Eer   = sqrt(mean((vEr - kron(ones(size(vEr,1),2),Er)).^2));
   Qe    = sqrt(mean((vQ - kron(ones(size(vQ,1),2),Q)).^2));

   Phimr = mean(vPhir);
   Phims = mean(vPhis);
   Emr   = mean(vEr);
   Ems   = mean(vEs);
   Qm    = mean(vQ);

   Phisr = std(vPhir);
   Phiss = std(vPhis);
   Esr   = std(vEr);
   Ess   = std(vEs);
   Qs    = std(vQ);

   mX1 = [Phims Ems Phimr Emr Qm];
   sX1 = [Phiss Ess Phisr Esr Qs];
   eX1 = [Phies Ees Phier Eer Qe];
