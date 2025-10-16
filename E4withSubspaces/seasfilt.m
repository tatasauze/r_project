function z1 = seasfilt(yt,s)

% Seasonal filter model with S(B)=1+B+B^2+...+B^{s-1}
% z = seasfilt(yt,s)

if any(isnan(yt))
    [th,di,la]=arma2thd([],[-1],[],[],1,s);
    yt = fismiss(th, di, yt);
end
S = ones(1,s-1);
[ts,ds,ls]=arma2thd([S],[],[],[],1,1);
z1 = residual(ts, ds, yt);      % Igual que fismod
z1 = z1(s:size(z1,1))/s; 