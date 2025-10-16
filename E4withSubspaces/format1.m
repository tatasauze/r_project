function  X = format1(vPhi,vE,vQ,vk)
%
   i = size(vk,2)/8;

   sQ = size(vQ,2)/i/8;
   sPhi = size(vPhi,2)/i/8;
   sE = size(vE,2)/i/8;

   Phim = reshape(mean(vPhi),8*sPhi,i)';
   Phis = reshape(std(vPhi),8*sPhi,i)';
   Em   = reshape(mean(vE),8*sE,i)';
   Es   = reshape(std(vE),8*sE,i)';
   Qm = reshape(mean(vQ),8*sQ,i)';
   Qs = reshape(std(vQ),8*sQ,i)';

   k = reshape([1:sQ^2],sQ,sQ);
   k = vech(k);

   X = zeros(i*2,sPhi*4+sE*8+(sQ+1)*sQ/2*8);

   for j=1:4
   %
       X(1:2:i*2,j:4:sPhi*4) = Phim(:,(j-1)*sPhi+1:j*sPhi);
       X(2:2:i*2,j:4:sPhi*4) = Phis(:,(j-1)*sPhi+1:j*sPhi);
   %

   k2 = sPhi*4;

   for j=1:8
   %
       X(1:2:i*2,(k2+j):8:(k2+sE*8)) = Em(:,(j-1)*sE+1:j*sE);
       X(2:2:i*2,(k2+j):8:(k2+sE*8)) = Es(:,(j-1)*sE+1:j*sE);
   %

   k2 = k2 + sE*8;

   for j=1:8
   %
       X(1:2:i*2,(k2+j):8:(k2+size(k,1)*8)) = Qm(:,(j-1)*sQ+k);
       X(2:2:i*2,(k2+j):8:(k2+size(k,1)*8)) = Qs(:,(j-1)*sQ+k);
   %