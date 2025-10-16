function [P, E, Q, P2] = Ricnit(Phi, H, G, L0)
%    [P, E, Q] = Ricnit(Phi, H, G, L0)
%    M�todo simpl�tico para la resoluci�n de la ec. de Riccatti

     n = size(Phi,1);
     iLam = inv(L0);

     A=[ Phi'-H'*iLam*G' zeros(n,n); -G*iLam*G'   eye(n) ];
     B=[eye(n) -H'*iLam*H; zeros(n,n) Phi-G*iLam*H ];

%     [W,d] = eig([ Phi'-H'*iLam*G' zeros(n,n);...
%                  -G*iLam*G'   eye(n) ], ...
%                   [eye(n) -H'*iLam*H; ...
%                   zeros(n,n) Phi-G*iLam*H ]);

     [W,d] = eig(A,B);

     d = diag(d);
     [e,index] = sort(abs(d));        % sort on magnitude of eigenvalues

%     if (~((e(n) < 1) & (e(n+1)>1)))
%           disp('Can''t order eigenvalues');
%     end

     WW = W(:,index);

     %%% State covariance:
     P = WW(n+1:2*n,1:n)*inv(WW(1:n,1:n));
     E = real((G - Phi * P * H')*inv(L0 - H * P * H'));
     Q = real(L0 - H * P * H');
     P2 = WW(n+1:2*n,n+1:2*n)*inv(WW(1:n,n+1:2*n));
