function [tchi2, desv, angle, rvchi, dg] = plotsay(y, i, s, ext)
%
%   Realiza un plot de los angulos principales
%
    if nargin < 4, ext = 0; end

    if nargin < 3, s = 1; end

    l = size(y,2);

    Yi = BlkHkl(y, 2*i, s, ext);
    Yj = Yi(i*l+1:(2*i)*l,:);
    Yi = Yi(1:i*l,:);
%    Yi = flipud(Yi(1:i*l,:));

    [Qi, R] = qr(Yi', 0);
    [Qj, R] = qr(Yj', 0);

    N = size(Qi,1);

    desv = zeros(l*i);
    tchi2 = zeros(l*i);

    for j = 1:i
    for k = 1:l
    %
        s = (j-1)*l+k;
        [U S V] = svd(Qi'*Qj(:,1:s));
        ang = diag(S(j:s,j:s));
        U = Qi * U(:,j:s);
        V = Qj(:,1:s) * V(:,j:s);

        sz = s-j+1;
        ds = zeros(sz,1);

        for h = 1:j-1
	%    
    	    ds = ds + diag((U(1:N-h,:)'*U(1+h:N,:)).*(V(1:N-h,:)'*V(1+h:N,:)));
	%
	end

        dv = 1 + 2*ds;
        desv(s,j:j+sz-1) = dv';
        vchi2 =  -N*flipud(cumsum(log(1-(ang(sz:-1:1).^2)./dv)));
        tchi2(s,j:j+sz-1) = chi2cdf(vchi2, (s-[1:sz]'-j+2).*(i*l-[1:sz]'-j+2))';
        rvchi(s,j:j+sz-1) = vchi2';
        dg(s,j:j+sz-1) = ((s-[1:sz]'-j+2).*(i*l-[1:sz]'-j+2))';
        angle(s,j:j+sz-1) = ang(1:sz)';
    %
    end
    end

    dv = sqrt(diag(desv))/sqrt(N);
    
    figure; whitebg('w'); close;
    figure;
    h1= axes('position',[0.1 0.55 0.8 0.35]);
    hold on;
    bar(diag(angle));
    plot(dv*2, '+');
    plot(dv*2, '-');
    set(h1, 'Box','on');
    title(['Angulos principales']);
  
    h2 = axes('position',[0.1 0.075 0.8 0.35]);
    hold on;
    bar(diag(tchi2));
    plot(ones(i*l,1)*.95, '+');
    set(h2, 'Box','on');
    title(['Probabilidad estadístico']);
    hold off
