function Z = BlkHkl(z, i, s, ext, inverse)

if nargin < 5, inverse = 0; end

if nargin < 4, ext = 0; end

if nargin < 3, s = 1; end

off = s*(i-1);

T = size(z,1);
l = size(z,2);

if ~ext
%
   Z = zeros(l*i, T-off);

   if ~inverse
       for j = 1:i
           Z(l*(j-1)+1:l*j,:) = z(s*(j-1)+1:T-s*(i-j),:)';
       end
   else
       for j = 1:i
           Z(l*(j-1)+1:l*j,:) = z(s*(i-j)+1:T-s*(j-1),:)';
       end
   end

%
else
%
   Z = zeros(l*i, T+off);

   if ~inverse
       for j = 1:i      
           Z(l*(j-1)+1:l*j,off+1:off+T) = z';
           off = off-s;
       end
   else
       for j = i:-1:1   
           Z(l*(j-1)+1:l*j,off+1:off+T) = z';
           off = off-s;
       end
   end
%
end