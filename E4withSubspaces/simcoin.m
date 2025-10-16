[t,d] = arma2thd([-2 1],[],[-.7],[],1,1);
N = 1000;
K = 20;
I = 6;
s1 = zeros(K,1);
s2 = zeros(K,1);

for i=1:K
i
    y = simmod(t,d,N);

    tst = diag(sidang(y,[],I));
    s1(i) = sum(tst >= .95);
    tst = diag(sidcoin(y,[],I));
    s2(i) = sum(tst >= .95);
end
