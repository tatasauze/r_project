function V = covcov(C,n,m,N)

V = zeros((n+1)*m*m);

for v = 1:m
    for s = 0:n
        for u = 1:m
            for j = v:m
                for h = 0:n
                    rw = (j-1)*(n+1)*m + h*m + 1;
                    cl = (v-1)*(n+1)*m + s*m + u; 
                    tsh = max(s - h, 0);
                    bsh = max(h - s, 0);
                    suma = zeros(m,1);
                    for t = -(n - bsh):(n - tsh)
                        at = abs(t);
                        x = (u-1)*(n+1)*m + at*m + 1;
                        at2 = abs(t+s-h);
                        y = (v-1)*(n+1)*m + at2*m + j; 
                        suma = suma + (1-at/N)*C(x:x+m-1)*C(y);
                    end
                    for t = -(n - h):(n - s)
                        at = abs(t);
                        at1 = abs(t+s); at2 = abs(t-h);
                        x = (v-1)*(n+1)*m + at1*m + 1;
                        y = (u-1)*(n+1)*m + at2*m + j; 
                        suma = suma + (1-at/N)*C(x:x+m-1)*C(y);
                    end
                    V(rw:rw+m-1,cl) = suma/N;
                end
            end
        end
    end
end

V2 = tril(V,-1);
V = diag(diag(V)) + V2 + V2';
