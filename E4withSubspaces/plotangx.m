function [ang, desv] = plotangx(y, u, i, s, ext)
%
%   Realiza un plot de los angulos principales
%
    if nargin < 5, ext = 0; end

    if nargin < 4, s = 1; end

    r = size(u,2);
    if r == 0
       [ang, desv] = plotang(y, i, s, ext);
       return;
    end

    l = size(y,2);

    Yi = [BlkHkl(u, 2*i, s, ext,1); BlkHkl(y, 2*i, s, ext)];
    Yj = Yi(i*(l+2*r)+1:i*2*(r+l),:);
    [Q, R] = qr(Yi(1:i*(l+2*r),:)', 0);
    Yi = R(r*i+1:(2*r+l)*i,r*i+1:(2*r+l)*i)'*Q(:,r*i+1:(2*r+l)*i)';
    
    index = zeros(i*(l+r),1);
    for k=1:i
    %
        index((k-1)*(r+l)+1:(k-1)*(r+l)+r) = (1:r)'+(k-1)*r;
        index((k-1)*(r+l)+r+1:k*(r+l)) = (1:l)'+(i-k)*l+i*r;
    %
    end

    Yi = Yi(index,:);

    [Qi, R] = qr(Yi', 0);
    [Qj, R] = qr(Yj', 0);

    N = min(size(Yi,2), size(y,1));

    desv = zeros(l*i,1);
    ang  = zeros(l*i,1);

    for j = 1:i
    %
        [U S V] = svd(Qi(:,1:j*(l+r))'*Qj(:,1:j*l),0);
        ang((j-1)*l+1:j*l) = diag(S((j-1)*l+1:j*l,(j-1)*l+1:j*l));
        U = Qi(:,1:j*(l+r)) * U(:,(j-1)*l+1:j*l);
        V = Qj(:,1:j*l) * V(:,(j-1)*l+1:j*l); 

        for k = 1:4
	%
     	    desv((j-1)*l+1:j*l,1) = desv((j-1)*l+1:j*l,1) + ...
             diag((U(1:N-k,:)'*U(1+k:N,:)).*(V(1:N-k,:)'*V(1+k:N,:)));
	%
	end
	
	desv((j-1)*l+1:j*l,1) = 2 * desv((j-1)*l+1:j*l,1) + 1;

    %
    end

    bar(ang);
    hold on;
    plot(sqrt(desv)*2/sqrt(size(y,1)), '+');
    plot(sqrt(desv)*2/sqrt(size(y,1)), '-');
    hold off;

%