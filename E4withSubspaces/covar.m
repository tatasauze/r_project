function Z = covar(x, e, l, XX)
%
%  x : Matriz (N*n) x (k*n)   (k regresores por ecuación)
%  e : Matriz n x N   (n ecuaciones)
%  l : número de retardos
%  la matriz de covarianzas es de dimensión k*n x k*n

   kn = size(x,2);
   n  = size(e,1);
   N  = size(e,2);

   X = zeros(kn*(l+1), 1);
   Y = zeros(kn*(l+1), kn);

N

   for t=1:N
       X(1:kn,:) = x((t-1)*n+1:t*n,:)'*e(:,t);
       Y = Y + X*X(1:kn,:)';
       X(kn+1:kn*(l+1),:) = X(1:kn*l,:);
   end

Y
   Z = Y(1:kn,:);

   for i=1:l
       Z = Z + Y(i*kn+1:(i+1)*kn,:) + Y(i*kn+1:(i+1)*kn,:)';
   end

Z
XX
   Z = XX*Z*XX;
