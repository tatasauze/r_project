function [ang, tchi2, desv] = plotang(y, i, s, ext)
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
    tchi2 = zeros(l*i,1);

    for j = 1:i
    %
        [U S V] = svd(Qi'*Qj(:,1:j*l));
        ang((j-1)*l+1:j*l) = diag(S((j-1)*l+1:j*l,(j-1)*l+1:j*l));
        U = Qi * U(:,(j-1)*l+1:j*l);
        V = Qj(:,1:j*l) * V(:,(j-1)*l+1:j*l);

        ds = zeros(l,1);

	k = 0;
	
        for k = 1:(j-1)
	%    
    	    ds = ds + diag((U(1:N-k,:)'*U(1+k:N,:)).*(V(1:N-k,:)'*V(1+k:N,:)));
	%
	end

        for h=1:l
	    if h > 1
               ds(h:l) = ds(h:l) + diag((U(1:N-k-h+1,h:l)'*U(1+k+h-1:N,h:l)).*(V(1:N-k-h+1,h:l)'*V(1+k+h-1:N,h:l)));
            end
  	    desv((j-1)*l+h:j*l) = 1 + 2*ds(h:l);
	    vchi2 =  -N*sum(log(1-(ang((j-1)*l+h:j*l).^2)./desv((j-1)*l+h:j*l)));
%            tchi2((j-1)*l+h) = chi2cdf(vchi2, (l+1-h)*(i*l-j+h));
            tchi2((j-1)*l+h) = chi2cdf(vchi2, (l+1-h)*((i-j+1)*l-h+1));
        end
    %
    end

    desv = sqrt(desv)/sqrt(N);
    
    figure; whitebg('w'); close;
    figure;
    h1= axes('position',[0.1 0.55 0.8 0.35]);
    hold on;
    bar(ang);
    plot(desv*2, '+');
    plot(desv*2, '-');
    set(h1, 'Box','on');
    title(['Angulos principales']);
   
    h2 = axes('position',[0.1 0.075 0.8 0.35]);
    hold on;
    bar(tchi2);
    plot(ones(i*l,1)*.95, '+');
    set(h2, 'Box','on');
    title(['Probabilidad estadístico']);
    hold off
%