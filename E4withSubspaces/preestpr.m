function theta2 = preestim(theta, din, z, userf)
%
% 7/3/97
% Copyright (c) Jaime Terceiro, 1997

global EEOPTION

if nargin < 3,  e4error(3); end
if nargin < 4,  userf = []; end

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ee(theta, din, userf);
N = size(z,1);
n = size(Phi,1);
m = size(H,1);
r = max([size(Gam,2), size(D,2)]);
if size(z,2) ~= m+r, e4error(11); end

if any(abs(eig(Phi)) >= EEOPTION(12)), stat = 0; else stat = 1; end

i = max(round(log(N)),ceil(n/m)+2);

if N < 2*i*(m+r+1-stat)-1, j = max(fix((N-i*(m+r+1)+1)/(m+r+1)),ceil(n/m)+1); else j = i; end
ij = i+j;

if N < ij*(m+r+1-stat)-1
   i = ceil(n/m)+1; j = i; ij = i+j;
   if N < i*(m+r+1-stat)+j*(r+1-stat)-1, error('No hay suficientes datos para utilizar preestim()'); end
end

if r
   Yi = blkhkel(z(:,1:m), ij, 1, stat);
   Ui = blkhkel(z(:,m+1:m+r), ij, 1, stat);
   N = min(size(Yi,2), size(z,1));
   [Q, R] = qr([Ui;Yi]', 0);
   R = R'/sqrt(N);
   ix = zeros(5,2);
   ix(:,1) = [1; j*r+1; ij*r+1; ij*r+i*m+1; ij*r+(i+1)*m+1];
   ix(:,2) = [ix(2:5,1)-1; ij*(m+r)];
else
   Yi = blkhkel(z, ij, 1, stat);
   N = min(size(Yi,2), size(z,1));
   [Q, R] = qr(Yi', 0);
   R = R'/sqrt(N);
   ix = zeros(3,2);
   ix(:,1) = [1; i*m+1; (i+1)*m+1];
   ix(:,2) = [ix(2:3,1)-1; ij*m];
end

if size(R,1) > size(R,2)
   R = [R zeros(size(R,1), size(R,1)-size(R,2))];
   pond = 0;
else
   pond = 1;
end

sth = size(theta,1);
if size(theta,2) == 1, theta2 = [theta zeros(sth,1)]; else, theta2 = theta; end

%if din(1,7) == 1 | din(1,7) == 3
%   theta2(sth-m+1:sth,2) = ones(m,1);
%else
%   theta2(sth-(m+1)*m/2+1:sth,2) = ones((m+1)*m/2,1);
%end

k = find(theta2(:,2) == 0);

if size(k,1) > 0
%   theta2(k,1) = zeros(size(k,1),1);

   %seteeopt('verbose','no');
   theta2 = eemin('preesti6', theta2, '', din, R, i+j, userf);
   sete4opt('verbose','si');
else
   pond = 0;
end

f = preesti6(theta2,din,R,i+j,userf);
%theta2 = eemin('preesti4', theta2, '', din, [K;B],theta2, [i j pond], userf);
return
if pond
   [f B] = preesti3(theta2,din,R,ix,[i j pond],userf);
else
   [f e] = fvarmax(theta2,din,z,userf);
   B = e'*e/size(z,1);
end

sth = size(theta2,1);

if din(1,7) == 1 | din(1,7) == 3
   if EEOPTION(5) == 1
      theta2(sth-m+1:sth,1) = sqrt(diag(B)); 
   else
      theta2(sth-m+1:sth,1) = diag(B);
   end
else
   if EEOPTION(5) == 1, B = cholp(B); end
   tmpndx = triu(ones(m));
   theta2(sth-(m+1)*m/2+1:sth,1) = B(tmpndx);
end

if size(theta,2) == 1, theta2 = theta2(:,1); else, theta2(:,2) = theta(:,2); end