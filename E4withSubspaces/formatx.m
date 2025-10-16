function  [mX, sX, eX]  = formatx(vPhi,vE,vQ,vG,vD,Phi,E,Q,D,G)
%
   i = size(vPhi,2)/size(Phi,1)/8;

   sQ = size(vQ,2)/i/8;
   sPhi = size(vPhi,2)/i/8;
   sE = size(vE,2)/i/8;
   sG = size(vG,2)/i/8;
   sD = size(vD,2)/i/8;


   Phi = Phi(:)';
   E   = E(:)';
   G   = G(:)';
   D   = D(:)';


   ssQ = size(Q,1);
   Q   = Q(:)';

   Phie = reshape(sqrt(mean((vPhi - kron(ones(size(vPhi,1),8*i),Phi)).^2)),8*sPhi,i)';
   Ee =   reshape(sqrt(mean((vE - kron(ones(size(vE,1),8*i),E)).^2)),8*sE,i)';
   Ge =   reshape(sqrt(mean((vG - kron(ones(size(vG,1),8*i),G)).^2)),8*sG,i)';
   De =   reshape(sqrt(mean((vD - kron(ones(size(vD,1),8*i),D)).^2)),8*sD,i)';
   Qe =   reshape(sqrt(mean((vQ - kron(ones(size(vQ,1),8*i),Q)).^2)),8*sQ,i)';
   Phim = reshape(mean(vPhi),8*sPhi,i)';
   Phis = reshape(std(vPhi),8*sPhi,i)';
   Em   = reshape(mean(vE),8*sE,i)';
   Es   = reshape(std(vE),8*sE,i)';
   Gm   = reshape(mean(vG),8*sG,i)';
   Gs   = reshape(std(vG),8*sG,i)';
   Dm   = reshape(mean(vD),8*sD,i)';
   Ds   = reshape(std(vD),8*sD,i)';
   Qm = reshape(mean(vQ),8*sQ,i)';
   Qs = reshape(std(vQ),8*sQ,i)';

   k = reshape([1:ssQ^2],ssQ,ssQ);
   k = vech(k);

%   mX = zeros(i,sPhi*4+sE*8+(ssQ+1)*ssQ/2*8);
   mX = zeros(i,sPhi*4+sG*4+sD*4+sE*8+sQ*8);
   sX = zeros(size(mX));
   eX = zeros(size(mX));

   for j=1:4
   %
       mX(:,j:4:sPhi*4) = Phim(:,(j-1)*sPhi+1:j*sPhi);
       sX(:,j:4:sPhi*4) = Phis(:,(j-1)*sPhi+1:j*sPhi);
       eX(:,j:4:sPhi*4) = Phie(:,(j-1)*sPhi+1:j*sPhi);

   %
   end

   k2 = sPhi*4;

   for j=1:4
   %
       mX(:,(k2+j):4:(k2+sD*4)) = Dm(:,(j-1)*sD+1:j*sD);
       sX(:,(k2+j):4:(k2+sD*4)) = Ds(:,(j-1)*sD+1:j*sD);
       eX(:,(k2+j):4:(k2+sD*4)) = De(:,(j-1)*sD+1:j*sD);
   %
   end

   k2 = k2 + sD*4;

   for j=1:4
   %
       mX(:,(k2+j):4:(k2+sG*4)) = Gm(:,(j-1)*sG+1:j*sG);
       sX(:,(k2+j):4:(k2+sG*4)) = Gs(:,(j-1)*sG+1:j*sG);
       eX(:,(k2+j):4:(k2+sG*4)) = Ge(:,(j-1)*sG+1:j*sG);
   %
   end

   k2 = k2 + sG*4;

   for j=1:4
   %
       mX(:,(k2+j):4:(k2+sE*4)) = Em(:,(j-1)*sE+1:j*sE);
       sX(:,(k2+j):4:(k2+sE*4)) = Es(:,(j-1)*sE+1:j*sE);
       eX(:,(k2+j):4:(k2+sE*4)) = Ee(:,(j-1)*sE+1:j*sE);
   %
   end

   k2 = k2 + sE*4;


   for j=1:4
   %
       mX(:,(k2+j):4:(k2+sQ*4)) = Qm(:,(j-1)*sQ+1:j*sQ);
       sX(:,(k2+j):4:(k2+sQ*4)) = Qs(:,(j-1)*sQ+1:j*sQ);
       eX(:,(k2+j):4:(k2+sQ*4)) = Qe(:,(j-1)*sQ+1:j*sQ);
   %
   end

   k2 = k2 + sQ*4;

   for j=1:4
   %
       mX(:,(k2+j):4:(k2+sE*4)) = Em(:,(j+3)*sE+1:(j+4)*sE);
       sX(:,(k2+j):4:(k2+sE*4)) = Es(:,(j+3)*sE+1:(j+4)*sE);
       eX(:,(k2+j):4:(k2+sE*4)) = Ee(:,(j+3)*sE+1:(j+4)*sE);
   %
   end

   k2 = k2 + sE*4;


   for j=1:4
   %
       mX(:,(k2+j):4:(k2+sQ*4)) = Qm(:,(j+3)*sQ+1:(j+4)*sQ);
       sX(:,(k2+j):4:(k2+sQ*4)) = Qs(:,(j+3)*sQ+1:(j+4)*sQ);
       eX(:,(k2+j):4:(k2+sQ*4)) = Qe(:,(j+3)*sQ+1:(j+4)*sQ);
   %
   end
