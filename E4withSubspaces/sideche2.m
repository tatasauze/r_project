function [Phi,sPhi,H,sH,E,sE,Qr,A1,A2,A3,A4,A5] = sidechel(y, u, i, ki, s)
%
% SIDECHEL -   Algoritmo de identificaci�n basado en subespacios, que permite
%             obtener una forma echelon forzando unos �ndices de observabilidad
%
% Modelos con variables ex�genas
%     [Phi,sPhi,H,sH,E,sE,Q,Gam,sGam,D,sD,innov] = sidechel(y, u, ij, ki, s)
%
% Modelos puramente estoc�sticos
%     [Phi,sPhi,H,sH,E,sE,Q,innov] = sidechel(y, [], ij, ki, s)
%
% 02/01/97


   global SIDOPTION

   if nargin < 5, s = 1; end

   exact = ~SIDOPTION(1);
   ext   = SIDOPTION(2);
   canon = SIDOPTION(3);
   pond  = SIDOPTION(4);
   oest  = SIDOPTION(5);

   if size(i,1) == 1, j = i(1,1); else, j = i(2,1); end
   i = i(1,1);

   n = sum(ki);

   if n > i | n > j, siderror(5); end

   m = size(y,2);
   mb = sum(ki > 0);
   r = size(u,2);

   if r
   %
     Yi = BlkHkel(y, i+j, s, ext);
     Ui = BlkHkel(u, i+j, s, ext);
     N = min(size(Yi,2), size(y,1));

% step 1
     [Q, R] = qr([Ui(r*(i+1)+1:r*(i+j),:);Ui(r*i+1:r*(i+1),:);Ui(1:r*i,:);Yi]', 0);
     R = R';
     ix = zeros(7,2);
     ix(:,1) = [1; (j-1)*r+1; j*r+1; (j+i)*r+1; (j+i)*r+m+1; (j+i)*r+i*m+1; (j+i)*r+(i+1)*m+1];
     ix(:,2) = [ix(2:7,1)-1; (i+j)*(m+r)];

% step 2
     if pond == 0
         Om = eye(j*m);
     elseif pond == 1
         Om = R(ix(6,1):ix(7,2),ix(6,1):ix(7,2))*R(ix(6,1):ix(7,2),ix(6,1):ix(7,2))'/N;
     else 
         Om = R(ix(6,1):ix(7,2),ix(3,1):ix(7,2))*R(ix(6,1):ix(7,2),ix(3,1):ix(7,2))';
     end

     sOm = sqrtm(Om);
     isOm = inv(sOm);

     if oest == 0
        if pond == 0
           Om1 = eye((j-1)*m);
        elseif pond == 1
           Om1 = R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))'/N;
        else 
           Om1 = R(ix(7,1):ix(7,2),ix(2,1):ix(7,2))*R(ix(7,1):ix(7,2),ix(2,1):ix(7,2))';
        end
        sOm1 = sqrtm(Om1);
        isOm1 = inv(sOm1);
     end

% step 3
     [U S V] = svd(isOm*R(ix(6,1):ix(7,2),ix(3,1):ix(5,2)));

     O = sOm * U(:,1:n);
     iO = U(:,1:n)'*isOm;

     if oest == 0
        [U1 S1 V1] = svd(isOm1*R(ix(7,1):ix(7,2),ix(2,1):ix(6,2)));
        O1 = sOm1 * U1(:,1:n);
        iO1 = pinv(O(1:(j-1)*m,:))*O1*U1(:,1:n)'*isOm1;
     else
        iO1  = pinv(O(1:(j-1)*m,:));
     end

% step 4
     Phi = iO1*O(m+1:j*m,:);
     H   = O(1:m,:);
     
% step 5
     [Phi, H, T] = echelon(ki, Phi, H);
     iT = inv(T);
     [strH, strPhi] = echelstr(ki);
     strH = strH(:);
     strPhi = strPhi(:);
     P = inv(R(ix(6,1):ix(6,2),ix(6,1):ix(6,2)));
     X = [iT*iO*R(ix(6,1):ix(7,2),ix(1,1):ix(5,2)) zeros(n,m)];
     Y = R(ix(6,1):ix(6,2),ix(1,1):ix(6,2));

     ik = find(ki > 0);
     YY = Y;
     YY(ik,:) = Y(ik,:) - X(1:mb,:);
     YY = P*YY;

     X2 = kron(X(1:mb,:)',P);
     idx = find(strH > 0);
     sidx = size(idx,1);
     X2 = [X2(:,idx) kron(R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))',P)];

     V = X2 \ YY(:);

     H1 = zeros(m*mb,1);
     H1(idx) = V(1:sidx);
     H1 = reshape(H1,m,mb);
     Di = reshape(V(sidx+1:size(V,1)),m,r*j);

     H = eye(m);
     H = H1 + H(:,ik);
     err = [YY(:)-X2*V];
     VV  = sqrt(diag((err'*err + m)/(m*N)*inv(X2'*X2)));
     H1(idx) = VV(1:sidx);
     sH = reshape(H1,m,mb);

     Zt = Y-H*X(1:mb,:)-Di*R(ix(1,1):ix(2,2),ix(1,1):ix(6,2));

     X2 = [kron(X',P) kron(R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))',P)];
     VV = (err'*err + m)/(m*N)*inv(X2'*X2)));
     VDi= VV(n*m+1:n*m*r*j,n*m+1:n*m*r*j);

     if n > mb
        Y = iT*iO1*R(ix(7,1):ix(7,2),ix(1,1):ix(6,2)) - Phi(:,mb+1:n)*X(mb+1:n,:);
     else
        Y = iT*iO1*R(ix(7,1):ix(7,2),ix(1,1):ix(6,2));
     end

     X2 = [X(1:mb,:) R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))];
     V = Y / X2;

     P = [Y-V*X2 iT*iO1*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))];
     P = inv(chol(P*P')');
     X2 = kron((H(ik,:)*X(1:mb,:))',P);

     YY = P*Y;
     idx = find(strPhi > 0);
     sidx = size(idx,1);

     X2 = [X2(:,idx) kron(R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))',P)];

     V = X2 \ YY(:);

     Phi1 = zeros(n*mb,1);
     Phi1(idx) = V(1:sidx);
     Phi(:,1:mb) = reshape(Phi1,n,mb);
     Gami = reshape(V(sidx+1:sidx+n*r*j), n, r*j);
     VV = inv(N*X2'*X2);
     VGami = VV(sidx+1:sidx+n*r*j,sidx+1:sidx+n*r*j);
     sPhi = zeros(n*mb,1);
     sPhi(idx) = sqrt(diag(VV(1:sidx,1:sidx)));
     sPhi = reshape(sPhi,n,mb);

     X2 = [H(ik,:)*X(1:mb,:);Zt+H*X(1:mb,:);R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))];

     V = Y / X2;
     E = V(:,mb+1:mb+m);
     E(:,ik) = -V(:,1:mb);
     err = [Y-V*X2 iT*iO1*R(ix(7,1):ix(7,2),ix(1,1):ix(6,2))];
     VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),mb+m+r*j,n)';
     sE = VV(:,mb+1:mb+m);
     sE(:,ik) = VV(:,1:mb);

     Qr =  Zt*Zt'/N;

     VGami = [VGami(n*j*(r-1)+1:n*j*r, n*j*(r-1)+1:n*j*r) VGami(n*j*(r-1)+1:n*j*r, 1:n*j*(r-1));...
              VGami(1:n*j*(r-1), n*j*(r-1)+1:n*j*r) VGami(1:n*j*(r-1), 1:n*j*(r-1))];
     VDi   = [VDi(m*j*(r-1)+1:m*j*r, m*j*(r-1)+1:m*j*r) VDi(m*j*(r-1)+1:m*j*r, 1:m*j*(r-1));...
              VDi(1:m*j*(r-1), m*j*(r-1)+1:m*j*r) VDi(1:m*j*(r-1), 1:m*j*(r-1))];

     VV = [VDi zeros(size(VDi)); zeros(size(VGami)) VGami];

     [A3 A1 A4 A2] = sidverg([[Di(:,j*(r-1)+1:j*r) Di(:,1:j*(r-1))] ; [Gami(:,j*(r-1)+1:j*r) Gami(:,1:j*(r-1))]], [eye(m) zeros(m,(j-1)*m); zeros(n,m) iT*iO1]*[eye(j*m)-O*iO], O, j, m, r, VV);
   
     if nargout > 6
        if ~ext
           A5 = Zt*Q(:,ix(1,1):ix(6,2))';
        else
           A5 = Zt*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(6,2))';
        end
     end
%
   else % r = 0, puramente estoc�sticos
   %
     Yi = BlkHkel(y, i+j, s, ext); 
     N = min(size(Yi,2), size(y,1));

% step 1
     [Q, R] = qr(Yi', 0);
     R = R';
     ix = zeros(3,2);
     ix(:,1) = [1; i*m+1; (i+1)*m+1];
     ix(:,2) = [ix(2:3,1)-1; (i+j)*m];

% step 2
     if pond == 0
         Om = eye(j*m);
     elseif pond == 1
         Om = R(ix(2,1):ix(3,2),ix(2,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(2,1):ix(3,2))'/N;
     else 
         Om = R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))';
     end

     sOm = sqrtm(Om);
     isOm = inv(sOm);

     if oest == 0
        if pond == 0
           Om1 = eye((j-1)*m);
        elseif pond == 1
           Om1 = R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))'/N;
        else 
           Om1 = R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))';
        end
        sOm1 = sqrtm(Om1);
        isOm1 = inv(sOm1);
     end

% step 3 & 4
     [U S V] = svd(isOm*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));

     O = sOm * U(:,1:n);
     iO = U(:,1:n)'*isOm;
     C = S(1:n,1:n) * V(:,1:n)'* R(ix(1,1):ix(1,2),ix(1,1):ix(1,2))'/ N;

     if oest == 0
        [U1 S1 V1] = svd(isOm1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2)));
        O1 = sOm1 * U1(:,1:n);
        iO1 = pinv(O(1:(j-1)*m,:))*O1*U1(:,1:n)'*isOm1;
     else
        iO1  = pinv(O(1:(j-1)*m,:));
     end

% step 5
     Phi = iO1*O(m+1:j*m,:);
     H   = O(1:m,:);

     [Phi, H, T] = echelon(ki, Phi, H);
     iT = inv(T);

%X = [iT*S(1:n,1:n)*V(:,1:n)' zeros(n,m)];
%Y = R(ix(2,1):ix(2,2),ix(1,1):ix(2,2));
%Zt = Y-H*X(1:mb,:);
%Y = iT*iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2)); 
%VV = Y / [X;Zt]
%res = [Y-VV*[X;Zt] iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))]*Q(:,ix(1,1):ix(3,2))';
%uidents(res',20);

     [strH, strPhi] = echelstr(ki);
     strH = strH(:);
     strPhi = strPhi(:);
     P = inv(R(ix(2,1):ix(2,2),ix(2,1):ix(2,2)));

     X = [iT*S(1:n,1:n)*V(:,1:n)' zeros(n,m)];
     Y = R(ix(2,1):ix(2,2),ix(1,1):ix(2,2));

     ik  = find(ki > 0);
     YY = Y;
     YY(ik,:) = Y(ik,:) - X(1:mb,:);
     YY = P*YY;

     X2 = kron(X(1:mb,:)',P);
     idx = find(strH > 0);
     X2 = X2(:,idx);

     V = X2 \ YY(:);

     H1 = zeros(m*mb,1);
     H1(idx) = V;
     H1 = reshape(H1,m,mb);
     H = eye(m);
     H = H1 + H(:,ik);
     err = [YY(:)-X2*V];
     VV  = sqrt(diag((err'*err + m)/(m*N)*inv(X2'*X2)));
%     VV  = sqrt(diag(inv(N*X2'*X2)));
     H1(idx) = VV;
     sH = reshape(H1,m,mb);

     Zt = Y-H*X(1:mb,:);

     if n > mb
        Y = iT*iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2)) - Phi(:,mb+1:n)*X(mb+1:n,:);
     else
        Y = iT*iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2));
     end

     P = [Y-Phi(:,1:mb)*X(1:mb,:) iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
     P = inv(chol(P*P')');
     X2 = kron((H(ik,:)*X(1:mb,:))',P);

     YY = P*Y;
     idx = find(strPhi > 0);
     sidx = size(idx,1);

     V = X2(:,idx) \ YY(:);

     Phi1 = zeros(n*mb,1);
     Phi1(idx) = V(1:sidx);
     Phi(:,1:mb) = reshape(Phi1,n,mb);
     sPhi = reshape(sqrt(diag(inv(N*X2'*X2))),n,mb);

     X2 = [H(ik,:)*X(1:mb,:);R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))];

     V = Y / X2;
     E = V(:,mb+1:mb+m);
     E(:,ik) = -V(:,1:mb);
     err = [Y-V*X2 iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
     VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),mb+m,n)';
     sE = VV(:,mb+1:mb+m);
     sE(:,ik) = VV(:,1:mb);
     
%err = [YY-P*Phi(:,1:mb)*X(1:mb,:)-P*E*Zt P*iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))]*Q(:,ix(1,1):ix(3,2))';
%Zv  = kron(Q(:,ix(1,1):ix(2,2))*Zt',P);
%Xt  = kron(Q(:,ix(1,1):ix(2,2))*X',P);

%(YY-P*Phi(:,1:mb)*X(1:mb,:)-P*E*Zt)*(YY-P*Phi*X-P*E*Zt)'/N
%P*iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))*(P*iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2)))'/N
%cov(err')
%Z = covar([Xt(:,idx) Zv], err, j-2, inv(X2'*X2))
%sqrt(diag(Z))

%     err = [YY(:)-X2*V];

%     VV  = sqrt(diag((err'*err + n)/(n*N)*inv(X2'*X2)));
%    VV  = sqrt(diag(inv(N*X2'*X2)));
%     Phi1(idx) = VV(1:sidx);
%     sPhi = reshape(Phi1,n,mb);
%     sE = reshape(VV(sidx+1:size(VV,1)),n,m);

     Qr =  Zt*Zt'/N;


% prueba
%     V = iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2)) / [T*X;Zt]
%     err = [iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2))-V*[T*X;Zt] iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
%     vPhi = kron(err*err'/N,S(1:n,1:n)^(-2))
%     vPhi = kron(T',iT)*vPhi*(kron(T',iT))'
%     sqrt(diag(vPhi))
% fin prueba          

    
     if nargout > 4
        if ~ext
           A1 = Zt*Q(:,ix(1,1):ix(2,2))';
        else
           A1 = Zt*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(2,2))';
        end
     end
   end
%
% end function




