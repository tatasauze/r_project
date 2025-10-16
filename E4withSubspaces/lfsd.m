function [f, z1, x0] = lfsd(theta, din, z)
% LFSD - Fast evaluation of the exact likelihood function for any time-invariant SS model
%    [f, innov, ssvect] = lfsd(theta, din, z)
% Under these conditions, this function outputs the same result that lfmod with
% vcond = 'idej'. It is valid for stationary and/or nonstationary models.
% theta  > parameter vector.
% din    > matrix which stores a description of the model dynamics.
% z      > matrix of observable variables.
% f      < value of the likelihood function.
% innov  < (optional) stores the sequence of innovations.
% ssvect < (optional) stores the sequence of values of the state vector.
%
% 7/6/07
% Copyright (c) Casals, Jerez and Sotoca, 2007

global E4OPTION

if nargin < 3,  e4error(3); end

saveinn = 0; if nargout > 1, saveinn = 1; end
scaleb  = E4OPTION(2);
econd   = E4OPTION(4);
zeps    = E4OPTION(15);

[H_D, type, m, r, s, n, np, userflag, userf, innov] = e4gthead(din);
[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
if econd == 5
   if r, econd = 2; else, econd = 3; end;
end

n = size(z,1);
l = size(Phi,1);
r = max([size(Gam,2), size(D,2)]);
if size(z,2) ~= m+r, e4error(11); end

if ~innov(1)
%
   [E,Q,U,iU] = sstoinn(Phi, E, H, C, Q, S, R);
%
else
%
   Q = C*Q*C';
   E = E * pinv(C);

   if scaleb, U = cholp(Q, abs(Q));
   else       U = cholp(Q); end
   iU = eye(m)/U';
%
end

if r & (econd == 1 | econd == 4)
%
   if econd == 1
      u0 = mean(z(:,m+1:m+r))';
   else
      u0 = z(1,m+1:m+r)';
   end
   [x0, Sigm, iSigm, nonstat, T,Z,RR,SS, uroots, ldet, RC, RS] = djccl(Phi, E*Q*E', 0,  Gam*u0);
%
else
   [x0, Sigm, iSigm, nonstat, T,Z,RR,SS, uroots, ldet, RC, RS] = djccl(Phi, E*Q*E', 0);   
end

WW  = zeros(l,l); WZ = zeros(l,1);
Phib = Phi - E*H;
Phibb0 = iU*H;

ff = 0.0;
uroots0 = uroots;

z1 = zeros(m,n);
Y = H;
YY = zeros(uroots);
T2 = T;

for t = 1:n  % main loop
    if r
    %
       z1(:,t)  = z(t,1:m)' - H*x0 - D*z(t,m+1:m+r)';
       x0  = Phi*x0 + Gam*z(t,m+1:m+r)' + E*z1(:,t);
    %
    else
    %
       z1(:,t)  = z(t,:)' - H*x0;
       x0  = Phi*x0 + E*z1(:,t);
    %
    end
    
    if uroots
       [HH, QZ, QT, rk, ign, k] = qr4(H*T2);
       YY = YY + T'*Y(k,:)'*Y(k,:)*T;
       Y = Y*Phi;
       T2 = Phi*T2*QT';
       uroots = uroots -rk;
    end
    WW  = WW + Phibb0'*Phibb0;
    WZ  = WZ + Phibb0'*iU*z1(:,t);
    Phibb0 = Phibb0*Phib;
%
end  % t

ff  =  2*n*sum(log(diag(U))) + sum(sum((iU*z1).^2));

if any(isnan([WW WZ])) | any(isinf([WW WZ])), e4error(25); end

WW = RS'*WW*RS;
WZ = RS'*WZ;

if ~nonstat  % stationary system
%
   Sigm = RC*Sigm*RC';
   if ~isempty(Sigm)
      M = cholp(Sigm);
      T = cholp(eye(size(M,1))+M*WW*M');
      ff = ff + 2*sum(log(diag(T)));
      if econd == 2
         T = cholp(WW);
         ff = ff - sum((T'\WZ).^2);
      else
         ff = ff - sum((T'\(M*WZ)).^2);
      end
   end
%
elseif nonstat == 2  % non stationary system
%
   [Ns S Ns] = svd(WW);
   S = diag(S);
   k = find(S > zeps);
   S = S(k);
   T = diag(1./sqrt(S))*Ns(:,k)';
   ff = ff + sum(log(S)) - sum((T*WZ).^2) - log(det(YY));
%
else  % partially stationary
%
   iSigm = RS'*iSigm*RS;
   [Ns S Ns] = svd(iSigm + WW);
   S = diag(S);
   k = find(S > zeps);
   S = S(k);
   T = diag(1./sqrt(S))*Ns(:,k)';
   ff = ff + sum(log(S)) - ldet - log(det(YY));
   if econd == 2
      [Ns S Ns] = svd(WW);    
      S = diag(S);
      k = find(S > zeps);
      S = S(k);
      T = diag(1./sqrt(S))*Ns(:,k)';
   end
   ff = ff - sum((T*WZ).^2);
%
end

f = 0.5*(ff + (n*m-uroots0)*log(2*pi));

if saveinn
%
   x0 = zeros(l,n+1);
   x0(:,1) = RS*pinv(WW)*WZ;
   for t = 1:n
       if r
       %
          z1(:,t)  = z(t,1:m)' - H*x0(:,t) - D*z(t,m+1:m+r)';
          x0(:,t+1)  = Phi*x0(:,t) + Gam*z(t,m+1:m+r)' + E*z1(:,t);
       %
       else
       %
          z1(:,t)  = z(t,:)' - H*x0(:,t);
          x0(:,t+1)  = Phi*x0(:,t) + E*z1(:,t);
       %
       end
   end

   z1 = z1';
   x0 = x0';
%
end
