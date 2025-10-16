function [Phi, Gam, E, H, D, C, Q, S, R] = mestaci(theta, din)

din = tomod(din);
np1 = din(1,6); m = din(1,2); r = din(1,3);

cf = theta(1,1);
theta = theta(2:size(theta,1),1);
nrows = m+1;
theta1 = theta(1:np1,1);  din1 = din(1:nrows,:);
[Phi, Gam, E, H, D, C, Q] = thd2ee(theta1, din1);
E(8,1) = cf;

nrows2 = size(din,1);
np2 = din(nrows+1,6); din2 = din(nrows+1:nrows2,:);
theta2 = theta(np1+1:size(theta,1),1);
[Phir, Gamr, Er, Hr, Dr, Cr, Qr] = thd2ee(theta2, din2);

l = size(Phi,1);
lf = size(Phir,1);
Phi = [Phi E*Hr; zeros(lf,l) Phir];
H   = [H Hr];
E   = [E; Er];
Q   = Qr;