function [f, innov, ssvect] = lfcd(theta, din, z)
% LFCD    - Computes the exact likelihood function of a general SS model.
%    [f, innov, ssvect] = lfcd(theta, din, z)
% theta    > parameter vector.
% din      > matrix which stores a description of the model dynamics.
% z        > matrix of observable variables.
% f        < value of the likelihood function.
% innov    < (optional) stores the sequence of innovations.
% ssvect   < (optional) stores the sequence of estimates of the state
%            vector.
%
% 7/6/07
% Copyright (c) Casals, Jerez and Sotoca, 2007


global E4OPTION

if nargin < 3,  e4error(3); end

saveinn = 0; if nargout >= 2, saveinn = 1; end
filtk   = 0; if E4OPTION(1) == 1, filtk = 1; end
scaleb  = E4OPTION(2);
vcond = E4OPTION(3);
econd = E4OPTION(4);
zeps  = E4OPTION(15);

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
n = size(z,1);
m = size(H,1);
r = size(Gam,2);

if econd == 5
   if r, econd = 2; else, econd = 3; end;
end

if econd == 2
   [f, innov, ssvect] = lffast(theta, din, z);
   return
end   

if size(z,2) ~= m+r, e4error(11); end

CRCt=C*R*C';  ESCt=E*S*C';  EQEt=E*Q*E';

if r & (econd == 1 | econd == 4)

   if econd == 1
      u0 = mean(z(:,m+1:m+r))';
   else
      u0 = z(1,m+1:m+r)';
   end

   [x0, P0, iSigm, nonstat, T,U2,RR,SS, uroots, ldet] = djccl(Phi, E*Q*E', 0,  Gam*u0);

else
   [x0, P0, iSigm, nonstat, T,U2,RR,SS, uroots, ldet] = djccl(Phi, E*Q*E', 0);
end


ff = 0.0;
uroots0 = uroots;

for t = 1:n  % main loop
%
    if r, z1  = z(t,1:m)' - H*x0 - D*z(t,m+1:m+r)';
    else  z1  = z(t,:)' - H*x0; end
    
    % Si existen raíces unitarias no identificadas se identifica la parte
    % de la innovación que no depende de esas raíces y que por tanto es
    % estacionaria
    if uroots
       % Descomposición QR con column pivoting
       [HH, QZ, QT, rk, ldetR12] = qr4(H*T);
    else
       rk = 0;
    end    
    % Filtro de Kalman estándar
    if saveinn, innov(t,:) = z1'; end
  
    B1  = H*P0*H' + CRCt;
    if scaleb, U = cholp(B1, abs(B1));
    else       U = cholp(B1); end
    iB1 = (eye(size(U))/U)/U';

    K1  = ((Phi*P0*H' + ESCt)/U)/U';
    if r, x1  = Phi*x0 + Gam*z(t,m+1:m+r)' + K1*z1;
    else  x1  = Phi*x0 + K1*z1; end
    Phib= Phi - K1*H;
    P1  = Phib*P0*Phib' + EQEt + K1*CRCt*K1' - K1*ESCt' - ESCt*K1';
    % Corrección al filtro de Kalman si ha raíces no estacionarias no
    % identificadas
    if uroots
       if rk
          WW = HH'*iB1*HH;
          if scaleb, U2 = cholp(WW, abs(WW));
          else       U2 = cholp(WW); end
          iWW = (eye(size(U2))/U2)/U2';
          WZ = HH'*iB1*z1;
          x1 = x1 + Phib*T*QZ'*iWW*WZ;          
          P1 = P1 + Phib*T*QZ'*iWW*QZ*T'*Phib';       
       end
       T = Phi*T*QT';
    end

    % A la función de verosimilitud sólo se añaden componentes
    % estacionarios
    if rk < m
       z1U = z1'/U;       
       ff  = ff + 2*sum(log(diag(U))) + z1U*z1U';
       if rk
          ff = ff - WZ'*iWW*WZ + 2*sum(log(diag(U2))) - ldetR12;
       end
    end
    
    % se restan las raíces unitarias identificadas en t
    uroots = uroots - rk;
    P0  = P1; x0  = x1; 
end  % t

f = 0.5*(ff + (n*m-uroots0)*log(2*pi));

if saveinn, [ff, innov, ssvect] = lffast(theta, din, z); end
