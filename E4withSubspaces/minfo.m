eeinit

[tm1, dm1,l1] = arma2thd([0],[],[.5],[],[1],1);
[tm2, dm2,l2] = arma2thd([0],[],[-.7],[],[1],1);
[tm3, dm3,l3] = arma2thd([0],[],[-.9],[],[1],1);

dst1 = [imod(tm1,dm1,zeros(100,1)) imod(tm1,dm1,zeros(200,1)) imod(tm1,dm1,zeros(1000,1))]
dst2 = [imod(tm2,dm2,zeros(100,1)) imod(tm2,dm2,zeros(200,1)) imod(tm2,dm2,zeros(1000,1))]
dst3 = [imod(tm3,dm3,zeros(100,1)) imod(tm3,dm3,zeros(200,1)) imod(tm3,dm3,zeros(1000,1))]

save imodma1 dst1 dst2 dst3

[tm1, dm1,l1] = arma2thd([-1 .64],[],[-.6 0],[],[1],1);
[tm2, dm2,l2] = arma2thd([0 .64],[],[-.6 0],[],[1],1);

dst1 = [imod(tm1,dm1,zeros(100,1)) imod(tm1,dm1,zeros(200,1)) imod(tm1,dm1,zeros(1000,1))]
dst2 = [imod(tm2,dm2,zeros(100,1)) imod(tm2,dm2,zeros(200,1)) imod(tm2,dm2,zeros(1000,1))]

save imodar2 dst1 dst2

[tm1, dm1] = arma2thd([0],[0],[-.5],[-.5],[1],12);
[tm2, dm2] = arma2thd([0],[0],[-.7],[-.9],[1],12);
[tm3, dm3] = arma2thd([0],[0],[-.9],[-.7],[1],12);

dst1 = [imod(tm1,dm1,zeros(100,1)) imod(tm1,dm1,zeros(200,1)) imod(tm1,dm1,zeros(1000,1))]
dst2 = [imod(tm2,dm2,zeros(100,1)) imod(tm2,dm2,zeros(200,1)) imod(tm2,dm2,zeros(1000,1))]
dst3 = [imod(tm3,dm3,zeros(100,1)) imod(tm3,dm3,zeros(200,1)) imod(tm3,dm3,zeros(1000,1))]

save imodsma1 dst1 dst2 dst3

[tm1, dm1] = arma2thd([-.5],[0],[0],[-.5],[1],12);
[tm2, dm2] = arma2thd([-.7],[0],[0],[-.9],[1],12);
[tm3, dm3] = arma2thd([-.9],[0],[0],[-.7],[1],12);

dst1 = [imod(tm1,dm1,zeros(100,1)) imod(tm1,dm1,zeros(200,1)) imod(tm1,dm1,zeros(1000,1))]
dst2 = [imod(tm2,dm2,zeros(100,1)) imod(tm2,dm2,zeros(200,1)) imod(tm2,dm2,zeros(1000,1))]
dst3 = [imod(tm3,dm3,zeros(100,1)) imod(tm3,dm3,zeros(200,1)) imod(tm3,dm3,zeros(1000,1))]

save imodsar1 dst1 dst2 dst3

[tm1, dm1] = arma2thd([-1 .64],[0],[0 0],[-.5],[1],12);
[tm2, dm2] = arma2thd([-1 .64],[0],[0 0],[-.7],[1],12);
[tm3, dm3] = arma2thd([-1 .64],[0],[0 0],[-.9],[1],12);

dst1 = [imod(tm1,dm1,zeros(100,1)) imod(tm1,dm1,zeros(200,1)) imod(tm1,dm1,zeros(1000,1))]
dst2 = [imod(tm2,dm2,zeros(100,1)) imod(tm2,dm2,zeros(200,1)) imod(tm2,dm2,zeros(1000,1))]
dst3 = [imod(tm3,dm3,zeros(100,1)) imod(tm3,dm3,zeros(200,1)) imod(tm3,dm3,zeros(1000,1))]

save imodsar2 dst1 dst2 dst3

[tm1, dm1,l1] = arma2thd([-.7 0  0;0 0 0;0 -.4 0],[],[0 1.1 0; 0 -.6 0; 0 0 .5],[],[1 -.7 .4;-.7 1 0;.4 0 1],1);
[tm2, dm2,l2] = arma2thd([-.4 -.3 .6;0 -.8 -.4;-.3 0 0],[],[-.7 0 0; -.1 -.2 0; .4 -.5 .1],[],[1 .5 .4;.5 1 .7;.4 .7 1],1);

dst1 = [imod(tm1,dm1,zeros(100,3)) imod(tm1,dm1,zeros(200,3)) imod(tm1,dm1,zeros(1000,3))]
dst2 = [imod(tm2,dm2,zeros(100,3)) imod(tm2,dm2,zeros(200,3)) imod(tm2,dm2,zeros(1000,3))]

save imodvma1 dst1 dst2

[te de] = arma2thd([],[],[],[],[.5],1);
[tm1, dm1] = arma2thd([-1 .64],[],[],[],[1],1);
[tm1,dm1] = comp2thd(tm1,dm1,[],[],te,de);
[tm2, dm2] = arma2thd([0 .64],[],[],[],[1],1);
[tm2,dm2] = comp2thd(tm2,dm2,[],[],te,de);

dst1 = [imod(tm1,dm1,zeros(100,1)) imod(tm1,dm1,zeros(200,1)) imod(tm1,dm1,zeros(1000,1))]
dst2 = [imod(tm2,dm2,zeros(100,1)) imod(tm2,dm2,zeros(200,1)) imod(tm2,dm2,zeros(1000,1))]

save imoderr dst1 dst2
