vPhi = zeros(200,2);
vE   = zeros(200,2);

for i=1:20
i
    y = simmod(t,d,50000);
i
    [vPhi(i,1),vPhi(i,2),H,sH,vE(i,1),vE(i,2),Q,innov] = sidechel(y, [], [3;3], 1);
end
    