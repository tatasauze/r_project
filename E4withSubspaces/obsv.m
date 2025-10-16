function  O = obsv(A,C)

n = size(A,1);
m = size(C,1);

O = zeros(m*n,n);
O(1:m,:) = C;
CA = C*A;

for i=2:n
    O((i-1)*m+1:i*m,:) = CA;
    CA = CA*A;
end
