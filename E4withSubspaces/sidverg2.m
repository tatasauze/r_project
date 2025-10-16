function [D, Gam, sD, sGam] = sidverg(Di, Gami, Op1, Op2, O, i, m, r, VOHd)
%
% SIDVERG -   Resuelve el sistema de ecuaciones propuesto por Viberg (1994)
%            para la identificación de las matrices D y Gam
%
%    [D, Gam] = sidverg(OHd, Op, O, i, m, r)
%
% 2/1/97


   n = size(O,2);
   l1 = size(Di,1);
   l2 = size(Gami,1);
   nm = n+m;

   YY1 = zeros(l*i,r);
   XX = zeros(l*i,nm);
   Ir = eye(r);
   X2 = zeros(r*l*i,r*nm);

   for k=1:i-1
       YY((k-1)*l+1:l*k,:) = OHd(:,(k-1)*r+1:r*k);
       XX((k-1)*l+1:l*k,:) = [Op(:,(k-1)*m+1:m*k) Op(:,m*k+1:m*i)*O(1:m*(i-k),:)];
       X2((k-1)*l*r+1:l*r*k,:) = kron(Ir,XX((k-1)*l+1:l*k,:));
   end

   YY((i-1)*l+1:l*i,:) = OHd(:,(i-1)*r+1:r*i);
   XX((i-1)*l+1:l*i,:) = [Op(:,(i-1)*m+1:m*i) zeros(l,n)];
   X2((i-1)*l*r+1:l*r*i,:) = kron(Ir,XX((i-1)*l+1:l*i,:));

   DGam = XX \ YY;
   D = DGam(1:m,:);
   Gam = DGam(m+1:m+n,:);

if nargin > 7
   X2 = pinv(X2);
   V = reshape(sqrt(diag(X2*VOHd*X2')),r,m+n)';
   sD = V(1:m,:);
   sGam = V(m+1:m+n,:);   
end