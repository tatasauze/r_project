sidinit
N = 200
sete4opt('econ','cero','vcon','cero');
sidopt('met','exac','obs','reestim','ext','si');


[t1,d1,l1] = tf2thd([],[],[-.7],[],1,1,[3; 2; 0; -4],[]);
[t2,d2,l2] = tf2thd([-1 .64],[],[-.6],[],1,1,[3; 2; 0; -4],[]);
[t3,d3,l3] = tf2thd([-.7],[],[],[],1,1,[3; 2; 0; -4],[]);

u = zeros(1000,4);
for i=1:4
  u(i:4:1000,i) = ones(250,1);
end

[M1, V1] = tabint(t1,d1,l1,u, [100;200;1000], N, [2:6], 1)
[M2, V2] = tabint(t2,d2,l2,u, [100;200;1000], N, [3:6], 2)
[M3, V3] = tabint(t3,d3,l3,u, [100;200;1000], N, [2:6], 1)
save inter1 M1 M2 M3 V1 V2 V3

