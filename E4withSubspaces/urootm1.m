function  [ur] = urootm1(N, i, S0); 

% Minimize the size given a minimum power (b): DEFAULT in NID
% corresponding to Gb(d > k)
% [ur] = urootm1(N, i, S0);

j = size(S0,1);

% Estimates of the penalty functions
C = [log(N)^2/N log(N)^2/N;
     exp(.43)*N^(-.39)*i^(-.07) exp(.43)*N^(-.39)*i^(-.07);
    -.3530+.03577*N-.00059*N^2+.0000030*N^3 exp(.18818)*N^(-.2854)*i^(-.17157);
    -.62117+.040609*N-.000556*N^2+.00000247*N^3 exp(1.5570)*N^(-.4694)*i^(-.4179);
    -.366001+.029814*N-.000301*N^2+.00000101*N^3-.063056*i exp(1.131294)*N^(-.360703)*i^(-.378167)];
 
Pos = zeros(j,1);  
if j>4
   for h = 1:4
       if N<88
       Pos(h,1) = 1-S0(h)^2-C(h,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(h,2);    
       end
   end
   for h = 5:j
       if N<121
       Pos(h,1) = 1-S0(h)^2-C(5,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(5,2);    
       end
   end
else 
   for h = 1:j
       if N<88
       Pos(h,1) = 1-S0(h)^2-C(h,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(h,2);    
       end
   end
end
 
ur = sum(Pos<0);