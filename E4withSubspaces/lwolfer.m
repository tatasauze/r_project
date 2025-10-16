%  LWOLFER.M - Ejemplo con la serie corta de Wolfer

!del wolferl.res
sidinit
diary wolferl.res
diary on
load wolferl.ser
%y = sqrt(wolferl);
y = wolferl;
y = y-mean(y);
sidopt('met','apro');

[tchi2, desv, angle, rvchi, dg] = sidang(y,[],4)
[tchi2, desv, angle, rvchi, dg] = sidang(y,[],6)
[Phi,sPhi,H,sH,E,sE,Q] = sidechei(y, [], 3:8, 2)
[Phi,H,E,Qr,innov,G,L0] = sidenti(y, [], 3:8, 2);
svd(obsv(Phi',E'))
[Phi, H, T, E] = echelon(2, Phi, H, E)
eig(Phi)
abs(eig(Phi))
G = inv(T)*G;
Q
[P, Q2, R] = siderv(Phi, H, Phi(:,1), G, L0)
P = lyapunov(Phi,Phi(:,1)*Q2*Phi(:,1)');
G = Phi*P*H'+Phi(:,1)*Q2;
L0 = H*P*H'+Q2+R;
[P, Et, Qt] = Ricnit(Phi, H, G, L0)

[tchi2, desv, angle, rvchi, dg] = sidang(innov',[],4)
[tchi2, desv, angle, rvchi, dg] = sidang(innov',[],6)

[Phi,sPhi,H,sH,E,sE,Q] = sidechei(y, [], 4:8, 3)
[Phi,H,E,Qr,innov,G,L0] = sidenti(y, [], 4:8, 3);
svd(obsv(Phi',E'))
[Phi, H, T, E] = echelon(3, Phi, H, E)
eig(Phi)
abs(eig(Phi))
G = inv(T)*G;
Q
[P, Q2, R] = siderv(Phi, H, Phi(:,1), G, L0)
P = lyapunov(Phi,Phi(:,1)*Q2*Phi(:,1)');
G = Phi*P*H'+Phi(:,1)*Q2;
L0 = H*P*H'+Q2+R;
[P, Et, Qt] = Ricnit(Phi, H, G, L0)


[tchi2, desv, angle, rvchi, dg] = sidang(innov',[],4)
[tchi2, desv, angle, rvchi, dg] = sidang(innov',[],6)
diary off