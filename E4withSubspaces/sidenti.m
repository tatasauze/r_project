function [Phix,Hx,Ex,Qx,A1x,A2x,A3x,A4x,A5x] = sidenti(y, u, i, n, s)
%
% Ver sident.m para sinopsis. La diferencia es que puede recibir en i un vector con
% posibles valores. Se seleccionará el valor de i que minimice la varianza del error de
% predicción

i = i(:);
r = size(u,2);

if nargin < 5, s = 1; end

fmax = inf;
for k=1:size(i,1)
   [Phi,H,E,Q,A1,A2,A3,A4,A5] = sident(y, u, i(k), n, s);
   if r
      [thet, di] = ee2thd(Phi, A1, E, H, A2, [], Q);
      f = pem(thet,di,[y u]);
   else
      [thet, di] = ee2thd(Phi, [], E, H, [], [], Q);
      f = pem(thet,di,y);
   end
   if f < fmax
      fmax = f;
      Phix = Phi;
      Hx   = H;
      Ex   = E;
      Qx   = Q;
      A1x = A1;
      A2x = A2;
      A3x = A3;
      A4x = A4;
      A5x = A5;
   end
end
