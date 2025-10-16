function [topt] = subseas(theta, din, z)
% SUBEST2   - Computes a fast estimate of the parameters of a model in THD format
%             with a subspace maximum likelihood function.
%           In innovations models it works with Covariance matrix 
%           In no-innovations models it works with Cholesky factors
%           [topt,std] = subest2(theta,din,z) % SUBEST2 CON ESTACIONALIDAD!!!
%
% theta  > parameter vector.
% din    > matrix which stores a description of the model dynamics.
% z      > matrix of observable variables.
% topt   < vector of estimates.
% std    < fast standard deviations of the estimates (inverse of the hessian or numerical approximation, at the optimum).
%
% March 2002
% Copyright (c) Jaime Terceiro, 1997

global E4OPTION

if nargin < 3,  e4error(3); end

[H_D, mtype, m, r, s, n, np, userflag, userf, innov, szpriv] = e4gthead(din);

if mtype < 100  % standard model
   if s~1; %% y arma model ;
       [FR, FS, AR, AS, V, G0] = thd2arm2(theta, din);
       [theta1, din1] = arma2thd([FS], [], [AS], [], 1, s);
       [toptS, hS] = subseas1(theta1, din1, z, 0);
       [FR_1, FS_1, AR_1, AS_1, V_1, G0_1] = thd2arm2(toptS, din1);       
       [theta1b, din1b] = arma2thd([], [FR_1], [], [AR_1], 1, s);
       z1 = residual(theta1b, din1b, z(:,1:m));
       [theta2, din2] = arma2thd([FR], [], [AR], [], 1, 1, G0, r);
       u = z(:,m+1:size(z,2));
       [topt, h] = subseas1(theta2, din2, [z1 u], 0);
       [fi, aa, ma, bb, var, mu] = thd2arm2(topt, din2);
       [topt, din] = arma2thd([fi], [FR_1], [ma], [AR_1], var, s, mu, r);
   else
       [topt, hess] = subes2(theta, din, z, 0);
   end
   %std = zeros(size(theta,1),1);  % pendiente de revision
   %if size(theta,2)>1
   %    std(find(theta(:,2)~=1)) = sqrt(diag(inv(hess)));
   %end
   return;
end

ptr1 = H_D+szpriv(2)+1;
[H_D, type, m, r, s, n, np, userflag, userf, innov, szpriv] = e4gthead(din(ptr1:size(din,1)));

if fix(mtype/100) == 1 % GARCH 
   theta1 = theta(1:np,:);
   din1 = din(ptr1:ptr1+H_D+szpriv(1)-1);
   topt   = subes1(theta1, din1, z(:,1:m+r), 0); % Modelo de la media 
   [f e]  = lffast(topt,din1,z(:,1:m+r));
   ptr1 = ptr1 + H_D + szpriv(1);
   din2   = din(ptr1:size(din,1));
   [H_D, type, mg, rg, s, n, npg, userflag, userf, innov, szpriv] = e4gthead(din2);
   theta2 = theta(np+1:size(theta,1),:);
   is_diag = rem(mtype,100);
   if is_diag
      ze = e.^2;
   else
      ze = zeros(size(e,1),mg);
      k = 1;
      for i=1:m
          for j=i:m
              ze(:,k) = e(:,i).*e(:,j);
              k=k+1;
          end
      end
   end
   if rg, ze = [ze z(:,m+r+1:m+r+rg)]; end
   topt = [topt; subes1(theta2, din2, ze, 1)]; % Modelo de la varianza
   topt = topt(1:size(theta,1),:);
   return;
end