function [Beta, VBeta] = simncov(Phi, Gam, E, H, D, C, Q, S, R, opt, T, N, p)
%
%  Simula el algoritmo JT para la estimación de matrices de covarianzas
%
%  [Beta, VBeta] = simncov(Phi, Gam, E, H, D, C, Q, S, R, opt, T, N)

[theta,din] = ss2thd(Phi, Gam, E, H, D, C, Q, S, R);

Beta = zeros(N,p);
VBeta= zeros(N,p);

sidopt('ext','no');

for j=1:N
%
j
    y = simmod(theta,din,T);
    [b, Vb] = sidvar(y, Phi, Gam, E, H, D, C, opt);
    Beta(j,:) = b';
    VBeta(j,:) = sqrt(diag(Vb))';
end
