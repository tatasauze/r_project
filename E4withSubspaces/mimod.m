function [Phi, Gam, E, H, D, C, Q, S, R] = mimod(theta, din)
[Phi, Gam, E, H, D, C, Q, S, R] = thd2ee(theta, tomod(din));