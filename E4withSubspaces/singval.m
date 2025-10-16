function [S] = singval(y, i, k, s)

% It computes canonical and no-canonical correlation (internal function) 
% S = singval(y, i, k, s)
% k = 1 -- No canonical
% k = 0 -- Canonical

if nargin < 4, s = 1; end
ext = 1; % Always with extended BHK matrices
j = i;   % dim(past) = dim(future)

[N m] = size(y);    
Yi = blkhkel(y, i+j, s, ext);
N = min(size(Yi,2), size(y,1));
     
[Q, R] = qr(Yi', 0);
R = R';
ix = zeros(3,2);
ix(:,1) = [1; i*m+1; (i+1)*m+1];
ix(:,2) = [ix(2:3,1)-1; (i+j)*m];

if k
  % No canonical
   [U S V] = svd(R(ix(2,1):ix(3,2),ix(2,1):ix(3,2)) \ R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));
else   
  % Canonical
   Om = R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))';
   isOm = pinv(sqrtm(Om));
   [U S V] = svd(isOm*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));       
end     
S = diag(S);