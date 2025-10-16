function [Phi, Gam, E, H, D, C, Q, S, R] = flores1(t,d)


d = tomod(d);
[Phi, Gam, E, H, D, C, Q, S, R] = thd2ee(t,d);
C = eye(3);
Gam=[];
D=[];