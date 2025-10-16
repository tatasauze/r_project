function [ang, desv] = plotang(y, i, s, ext)
%
%   Realiza un plot de los angulos principales
%
    if nargin < 4, ext = 0; end

    if nargin < 3, s = 1; end

    l = size(y,2);

    Yi = BlkHkl(y, 2*i, s, ext);
    Yj = Yi(i*l+1:(2*i)*l,:);
    Yi = flipud(Yi(1:i*l,:));

    [Qi, R] = qr(Yi', 0);
    [Qj, R] = qr(Yj', 0);

    N = size(Qi,1);

    desv = zeros(l*i,1);
    ang  = zeros(l*i,1);

    for j = 1:i
    %
        [U S V] = svd(Qi(:,1:j*l)'*Qj(:,1:j*l));        
        ang((j-1)*l+1:j*l) = diag(S((j-1)*l+1:j*l,(j-1)*l+1:j*l));
        U = Qi(:,1:j*l) * U(:,(j-1)*l+1:j*l);
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

    desv = sqrt(desv)*2/sqrt(size(y,1));

%    bar(ang);
%    hold on;
%    plot(sqrt(desv)*2/sqrt(size(y,1)), '+');
%    plot(sqrt(desv)*2/sqrt(size(y,1)), '-');
%    hold off;

%