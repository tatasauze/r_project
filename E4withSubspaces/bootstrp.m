function [vPhi, vE, vQ] = Bootstrp(y, i, n, s)


if nargin < 4,   s = 1; end

N = 100;

[Phi,H,E,Q,innov] = sident(y, [], i, n, s);
[Phiar, Har, Tr, Ear] = echelon(n,Phi, H, E)
n = size(Phi,1);
m = size(H,1);
Ear = Ear + Phiar(:,1:m);

T = size(y,1);
Ti = size(innov,2)-1;
rand('seed',sum(100*clock));

vPhi = zeros(N,n*m);
vE   = zeros(N,n*m);
vQ   = zeros(N,m*m);

for j=1:N
%
    k = round(rand(T,1)*Ti + 1);
    a = innov(:,k);
    x = zeros(n,1);

    for h=1:T
    %
         y(h,:) = (Har*x + a(:,h))';
	 x = Phiar * x + Ear*a(:,h);
    %
    end

    [Phi,H,E,Q] = sident(y, [], i, n, s);
    [Phi, H, Tr, E] = echelon(n,Phi, H, E);
    F = Phi(:,1:m);

    vPhi(j,:) = F(:)';
    vE(j,:)   = E(:)';
    vQ(j,:)   = Q(:)';
%
end
    