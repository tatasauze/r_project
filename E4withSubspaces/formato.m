function  [mX, sX, eX]  = formato(vPhi,vE,vQ,vk,Phi,E,Q)
%
   i = size(vk,2)/8;

   sQ = size(vQ,2)/i/8;
   sPhi = size(vPhi,2)/i/8;
   sE = size(vE,2)/i/8;

   Phi = Phi(:)';
   E   = E(:)';

   ssQ = size(Q,1);
   Q   = Q(:)';

   Phie = reshape(sqrt(mean((vPhi - kron(ones(size(vPhi,1),8*i),Phi)).^2)),8*sPhi,i)';
   Ee =   reshape(sqrt(mean((vE - kron(ones(size(vE,1),8*i),E)).^2)),8*sE,i)';
   Qe =   reshape(sqrt(mean((vQ - kron(ones(size(vQ,1),8*i),Q)).^2)),8*sQ,i)';
   Phim = reshape(mean(vPhi),8*sPhi,i)';
    Phis = reshape(std(vPhi),8*sPhi,i)';
   Em   = reshape(mean(vE),8*sE,i)';
   Es   = reshape(std(vE),8*sE,i)';
   Qm = reshape(mean(vQ),8*sQ,i)';
   Qs = reshape(std(vQ),8*sQ,i)';

   k = reshape([1:ssQ^2],ssQ,ssQ)
   k = vech(k)

%   mX = zeros(i,sPhi*4+sE*8+(ssQ+1)*ssQ/2*8);
   mX = zeros(i,sPhi*4+sE*8+sQ*8);
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

   for j=1:8
   %
       mX(:,(k2+j):8:(k2+sE*8)) = Em(:,(j-1)*sE+1:j*sE);
       sX(:,(k2+j):8:(k2+sE*8)) = Es(:,(j-1)*sE+1:j*sE);
       eX(:,(k2+j):8:(k2+sE*8)) = Ee(:,(j-1)*sE+1:j*sE);
   %
   end

   k2 = k2 + sE*8;

   for j=1:8
   %
%       mX(:,(k2+j):8:(k2+size(k,1)*8)) = Qm(:,(j-1)*sQ+k);
%       sX(:,(k2+j):8:(k2+size(k,1)*8)) = Qs(:,(j-1)*sQ+k);
%       eX(:,(k2+j):8:(k2+size(k,1)*8)) = Qe(:,(j-1)*sQ+k);
       mX(:,(k2+j):8:(k2+sQ*8)) = Qm(:,(j-1)*sQ+1:j*sQ);
       sX(:,(k2+j):8:(k2+sQ*8)) = Qs(:,(j-1)*sQ+1:j*sQ);
       eX(:,(k2+j):8:(k2+sQ*8)) = Qe(:,(j-1)*sQ+1:j*sQ);
   %
   end