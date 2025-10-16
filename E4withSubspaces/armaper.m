function [Phi, H, T, E] = armaper(k, A, C, E)
%
   s = size(k,1);
   m = size(C,1) / s;
   sk = sum(k);
   T = zeros(sk, max(k));
   k = [k;k(1)];
   pos = cumsum([1;k(2:s+1)]);
   post = cumsum([1;k(1:s)]);

   Phi = zeros(size(A));
   H = zeros(size(C));

   Nk1 = k(s);
   for i=1:s
   %   
       nt = k(i) - min(k(i)+m, Nk1) + m
       if nt > 0
          O = zeros(k(i));
          M1 = eye(k(i));
          for j=1:k(i)/m
              h = rem(i+j-1,s);
              if h == 0, h = s; end
	      O((j-1)*m+1:j*m,:) = C((h-1)*m+1:h*m, 1:k(h)) * M1;
              M1 = A(pos(h):pos(h+1)-1, 1:k(h))*M1;
          end
O
          iO = inv(O)
          T(post(i):post(i+1)-1,k(i)-nt+1:k(i)) = iO(:,k(i)-nt+1:k(i))
       end
       if i < s 
          j = min(k(i+1)+m,k(i));
          if j > m
             T(post(i+1):post(i+2)-1,1:j-m) = A(pos(i):pos(i+1)-1, 1:k(i))*T(post(i):post(i+1)-1,m+1:j);
          end
       else
          j = min(k(1)+m,k(i));
          if j > m
             T(1:k(1),1:j-m) = A(pos(i):pos(i+1)-1, 1:k(i))*T(post(s):post(s+1)-1,m+1:j);
          end
       end
       Nk1 = k(i);
   end

   for i=1:s-1
   %   
       if i < s 
          j = min(k(i+1)+m,k(i));
          if j > m
             T(post(i+1):post(i+2)-1,1:j-m) = A(pos(i):pos(i+1)-1, 1:k(i))*T(post(i):post(i+1)-1,m+1:j);
          end
       else
          j = min(k(1)+m,k(i));
          if j > m
             T(1:k(1),1:j-m) = A(pos(i):pos(i+1)-1, 1:k(i))*T(post(s):post(s+1)-1,m+1:j);
          end
       end
   end
T

   for i=1:s
   %   
       if i == s
          iT = inv(T(1:k(1),1:k(1)));
       else
          iT = inv(T(post(i+1):post(i+2)-1,1:k(i+1)));
       end
       Phi(pos(i):pos(i+1)-1,1:k(i)) = iT * A(pos(i):pos(i+1)-1,1:k(i)) * T(post(i):post(i+1)-1,1:k(i));
       E(pos(i):pos(i+1)-1,:) = iT * E(pos(i):pos(i+1)-1,:) - Phi(pos(i):pos(i+1)-1,1:m);
       H((i-1)*m+1:i*m, 1:k(i)) = C((i-1)*m+1:i*m, 1:k(i)) * T(post(i):post(i+1)-1,1:k(i));
   %
   end
