function  [ur] = urootm(N, i, S0); 

% Minimize the size with a consisten power (a)
% corresponding to Ga(d > k) 
% [ur] = urootm(N, i, S0);

j = size(S0,1);

% Estimates of the penalty functions
C = [exp(.61)*N^(-.495)*i^(-.10) exp(.61)*N^(-.495)*i^(-.10);
     exp(.67)*N^(-.39)*i^(-.06) exp(.67)*N^(-.39)*i^(-.06);
    -.30532+.03968*N-.00065*N^2+.0000033*N^3 exp(.78648)*N^(-.32804)*i^(-.2259); 
    -.634739+.044421*N-.000603*N^2+.00000267*N^3 exp(1.589335)*N^(-.436601)*i^(-.364576);
    -.316565+.031558*N-.000315*N^2+.00000105*N^3-.075574*i exp(1.312803)*N^(-.38284)*i^(-.280176)];

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