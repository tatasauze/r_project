function [dPhi,dGam,dE,dH,dD,dC,dQ,dS,dR] = gflore1(t, d, i)

d = tomod(d);
[dPhi,dGam,dE,dH,dD,dC,dQ,dS,dR] = ee_dv(t, d, i);
dC = zeros(3);
dGam=[];
dD=[];